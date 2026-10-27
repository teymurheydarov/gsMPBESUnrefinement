/** @file gsMPBES.hpp

    @brief Implementation of gsMPBES

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): T. Heydarov
*/

#pragma once

#include <gsMPBESUnrefinement/src/gsMPBES.h>
#include <gsCore/gsAffineFunction.h>
#include <gsUtils/gsPointGrid.h>

#include <algorithm>
#include <map>
#include <numeric>
#include <random>
#include <vector>

namespace gismo
{

template<class T>
gsMPBES<T>::gsMPBES(const gsMultiBasis<T> & thb)
: Base(), m_numLevels(0)
{
    init(thb);
}

template<class T>
gsMPBES<T>::gsMPBES(const gsMultiPatch<T> & mp)
: Base(), m_numLevels(0)
{
    if ( mp.nInterfaces() == 0 && mp.nBoundary() == 0 )
    {
        gsMultiPatch<T> copy(mp);
        copy.computeTopology();
        init( gsMultiBasis<T>(copy) );
    }
    else
        init( gsMultiBasis<T>(mp) );
}

template<class T>
const gsTHBSplineBasis<2,T> & gsMPBES<T>::thbBasis(index_t p) const
{
    const gsTHBSplineBasis<2,T> * b = dynamic_cast<const gsTHBSplineBasis<2,T>*>(&m_thb.basis(p));
    GISMO_ENSURE(b, "gsMPBES: patch "<<p<<" does not carry a gsTHBSplineBasis<2>.");
    return *b;
}

template<class T>
void gsMPBES<T>::init(const gsMultiBasis<T> & thb)
{
    GISMO_ENSURE(thb.nBases() != 0, "gsMPBES: no patches.");
    m_thb = thb;
    for (size_t p = 0; p != m_thb.nBases(); ++p)
        thbBasis(p); // type check
    buildLevels();
    identifyTwins();
    selectFunctions();
    truncateAndMap();
}

template<class T>
void gsMPBES<T>::buildLevels()
{
    const size_t nP = m_thb.nBases();
    m_levels.resize(nP);
    m_numLevels = 0;
    for (size_t p = 0; p != nP; ++p)
    {
        const gsTHBSplineBasis<2,T> & thb = thbBasis(p);
        const index_t nL = static_cast<index_t>(thb.maxLevel()) + 1;
        m_numLevels = math::max(m_numLevels, nL);
        // the levels are the dyadic refinements of level 0 (as the tree of
        // the THB basis assumes)
        gsTensorBSplineBasis<2,T> b = thb.tensorLevel(0);
        m_levels[p].clear();
        for (index_t l = 0; l != nL; ++l)
        {
            m_levels[p].push_back(b);
            if (l + 1 != nL) b.uniformRefine();
        }
        GISMO_ENSURE(m_levels[p][nL-1].size() == thb.tensorLevel(nL-1).size(),
                     "gsMPBES: the levels of patch "<<p<<" are not dyadic refinements of level 0.");
    }
    // the same degree everywhere
    const short_t d0 = m_levels[0][0].degree(0), d1 = m_levels[0][0].degree(1);
    for (size_t p = 0; p != nP; ++p)
        GISMO_ENSURE(m_levels[p][0].degree(0) == d0 && m_levels[p][0].degree(1) == d1,
                     "gsMPBES: patch "<<p<<" has a different degree.");

    // flat index spaces of the levels
    m_levelOffset.assign(m_numLevels, std::vector<index_t>(nP, 0));
    m_patchOfFlat.assign(m_numLevels, std::vector<index_t>());
    m_globalIndex.resize(nP);
    for (index_t l = 0; l != m_numLevels; ++l)
    {
        index_t off = 0;
        for (size_t p = 0; p != nP; ++p)
        {
            m_levelOffset[l][p] = off;
            if (l <= maxLevel(p))
            {
                off += m_levels[p][l].size();
                m_patchOfFlat[l].insert(m_patchOfFlat[l].end(), m_levels[p][l].size(), static_cast<index_t>(p));
            }
        }
    }
    for (size_t p = 0; p != nP; ++p)
    {
        m_globalIndex[p].resize(m_levels[p].size());
        for (size_t l = 0; l != m_levels[p].size(); ++l)
            m_globalIndex[p][l].assign(m_levels[p][l].size(), -1);
    }
}

template<class T>
void gsMPBES<T>::identifyTwins()
{
    const size_t nP = m_thb.nBases();
    const gsBoxTopology & topol = m_thb.topology();
    m_classId.resize(m_numLevels);
    m_classes.resize(m_numLevels);
    m_noTwin.resize(m_numLevels);
    m_insideMemo.resize(m_numLevels);

    for (index_t l = 0; l != m_numLevels; ++l)
    {
        index_t total = 0;
        for (size_t p = 0; p != nP; ++p)
            if (l <= maxLevel(p)) total += m_levels[p][l].size();
        std::vector<index_t> parent(total);
        std::iota(parent.begin(), parent.end(), 0);
        m_noTwin[l].assign(total, 0);
        m_insideMemo[l].assign(total, -1);

        auto find = [&](index_t i)
        {
            while (parent[i] != i) { parent[i] = parent[parent[i]]; i = parent[i]; }
            return i;
        };
        auto unite = [&](index_t a, index_t b)
        {
            a = find(a); b = find(b);
            if (a != b) parent[math::max(a,b)] = math::min(a,b);
        };

        for (gsBoxTopology::const_iiterator it = topol.iBegin(); it != topol.iEnd(); ++it)
        {
            const boundaryInterface & bi = *it;
            const index_t p1 = bi.first().patch, p2 = bi.second().patch;
            const bool has1 = (l <= maxLevel(p1)), has2 = (l <= maxLevel(p2));
            if (has1 && has2)
            {
                const gsTensorBSplineBasis<2,T> & b1 = m_levels[p1][l];
                const gsTensorBSplineBasis<2,T> & b2 = m_levels[p2][l];
                gsMatrix<index_t> bnd1, bnd2;
                bnd1 = b1.boundary(bi.first().side());
                bnd2 = b2.boundary(bi.second().side());
                GISMO_ENSURE(bnd1.rows() == bnd2.rows(),
                             "gsMPBES: the knot vectors of patches "<<p1<<" and "<<p2
                             <<" do not match along their interface on level "<<l<<".");
                b1.matchWith(bi, b2, bnd1, bnd2);
                for (index_t i = 0; i != bnd1.rows(); ++i)
                    unite( flat(p1, l, bnd1(i,0)), flat(p2, l, bnd2(i,0)) );
            }
            else if (has1 || has2)
            {
                // the neighbour has no level l: the boundary functions of
                // this side cannot be part of a C0 function of level l
                const index_t p = has1 ? p1 : p2;
                const boxSide s = has1 ? bi.first().side() : bi.second().side();
                const gsMatrix<index_t> bnd = m_levels[p][l].boundary(s);
                for (index_t i = 0; i != bnd.rows(); ++i)
                    m_noTwin[l][ flat(p, l, bnd(i,0)) ] = 1;
            }
        }

        // classes from the roots
        m_classId[l].assign(total, -1);
        m_classes[l].clear();
        std::map<index_t, index_t> idOfRoot;
        for (index_t i = 0; i != total; ++i)
        {
            const index_t r = find(i);
            typename std::map<index_t, index_t>::iterator it = idOfRoot.find(r);
            if (it == idOfRoot.end())
            {
                it = idOfRoot.insert(std::make_pair(r, static_cast<index_t>(m_classes[l].size()))).first;
                m_classes[l].push_back(std::vector<index_t>());
            }
            m_classId[l][i] = it->second;
            m_classes[l][it->second].push_back(i);
        }
    }
}

template<class T>
int gsMPBES<T>::componentLevel(index_t p, index_t l, index_t index) const
{
    if (l > maxLevel(p)) return -1;
    if (m_noTwin[l][ flat(p, l, index) ]) return -1;
    const gsTensorBSplineBasis<2,T> & b = m_levels[p][l];
    const gsVector<index_t,2> ti = b.tensorIndex(index);
    gsVector<index_t,2> low, upp;
    gsMatrix<index_t,1,2> es;
    for (short_t k = 0; k != 2; ++k)
    {
        b.component(k).elementSupport_into(ti[k], es);
        low[k] = es(0,0);
        upp[k] = es(0,1);
    }
    return thbBasis(p).tree().query3(low, upp, static_cast<int>(l));
}

template<class T>
bool gsMPBES<T>::classInside(index_t p, index_t l, index_t index) const
{
    const index_t f = flat(p, l, index);
    signed char & memo = m_insideMemo[l][f];
    if (memo >= 0) return memo != 0;
    bool inside = true;
    const std::vector<index_t> & cls = m_classes[l][ m_classId[l][f] ];
    for (size_t i = 0; i != cls.size() && inside; ++i)
    {
        const index_t q   = m_patchOfFlat[l][cls[i]];
        const index_t idx = cls[i] - m_levelOffset[l][q];
        if ( componentLevel(q, l, idx) < static_cast<int>(l) ) inside = false;
    }
    memo = inside ? 1 : 0;
    return inside;
}

template<class T>
void gsMPBES<T>::selectFunctions()
{
    m_components.clear();
    m_level.clear();
    for (index_t l = 0; l != m_numLevels; ++l)
    {
        for (size_t c = 0; c != m_classes[l].size(); ++c)
        {
            const std::vector<index_t> & cls = m_classes[l][c];
            bool allInside = true, someExact = false;
            std::vector<Component> comps;
            for (size_t i = 0; i != cls.size(); ++i)
            {
                Component cp;
                cp.patch = m_patchOfFlat[l][cls[i]];
                cp.level = l;
                cp.index = cls[i] - m_levelOffset[l][cp.patch];
                const int m = componentLevel(cp.patch, l, cp.index);
                if (m < static_cast<int>(l)) allInside = false;
                if (m == static_cast<int>(l)) someExact = true;
                comps.push_back(cp);
            }
            if ( !(allInside && someExact) ) continue;
            const index_t g = static_cast<index_t>(m_components.size());
            for (size_t i = 0; i != comps.size(); ++i)
                m_globalIndex[comps[i].patch][l][comps[i].index] = g;
            m_components.push_back(comps);
            m_level.push_back(l);
        }
    }
    m_truncated.assign(m_components.size(), false);
}

template<class T>
void gsMPBES<T>::truncateAndMap()
{
    const size_t nP = m_thb.nBases();

    // two-scale matrices between consecutive levels of every patch
    std::vector< std::vector< gsSparseMatrix<T,RowMajor> > > transfer(nP);
    for (size_t p = 0; p != nP; ++p)
    {
        transfer[p].resize(m_levels[p].size());
        for (size_t k = 1; k != m_levels[p].size(); ++k)
        {
            gsTensorBSplineBasis<2,T> b = m_levels[p][k-1];
            b.uniformRefine_withTransfer(transfer[p][k], 1, 1);
            GISMO_ENSURE(b.size() == m_levels[p][k].size() && transfer[p][k].rows() == m_levels[p][k].size(),
                         "gsMPBES: inconsistent refinement of level "<<k<<" on patch "<<p<<".");
        }
    }

    // offsets of the finest-level bases in the local index space
    std::vector<index_t> fineOffset(nP + 1, 0);
    for (size_t p = 0; p != nP; ++p)
        fineOffset[p+1] = fineOffset[p] + m_levels[p][maxLevel(p)].size();

    std::vector< gsEigen::Triplet<T> > triplets;
    gsMatrix<T> v, w;
    for (size_t g = 0; g != m_components.size(); ++g)
    {
        for (size_t ci = 0; ci != m_components[g].size(); ++ci)
        {
            const Component & c = m_components[g][ci];
            const index_t Kp = maxLevel(c.patch);
            v.setZero(m_levels[c.patch][c.level].size(), 1);
            v(c.index, 0) = 1;
            for (index_t k = c.level + 1; k <= Kp; ++k)
            {
                w.noalias() = transfer[c.patch][k] * v;
                for (index_t j = 0; j != w.rows(); ++j)
                    if ( w(j,0) != 0 && classInside(c.patch, k, j) )
                    {
                        w(j,0) = 0;
                        m_truncated[g] = true;
                    }
                v.swap(w);
            }
            for (index_t j = 0; j != v.rows(); ++j)
                if ( v(j,0) != 0 )
                    triplets.push_back( gsEigen::Triplet<T>(fineOffset[c.patch] + j, static_cast<index_t>(g), v(j,0)) );
        }
    }
    gsSparseMatrix<T> M(fineOffset[nP], static_cast<index_t>(m_components.size()));
    M.setFromTriplets(triplets.begin(), triplets.end());

    std::vector< gsBasis<T>* > fine;
    for (size_t p = 0; p != nP; ++p)
        fine.push_back( new gsTensorBSplineBasis<2,T>(m_levels[p][maxLevel(p)]) );
    gsMultiBasis<T> fineBases(fine, m_thb.topology());
    Base::init(fineBases, M);
}

template<class T>
index_t gsMPBES<T>::globalIndex(index_t p, index_t l, index_t index) const
{
    if (l > maxLevel(p)) return -1;
    return m_globalIndex[p][l][index];
}

template<class T>
index_t gsMPBES<T>::numFunctionsOfLevel(index_t l) const
{
    return static_cast<index_t>(std::count(m_level.begin(), m_level.end(), l));
}

template<class T>
index_t gsMPBES<T>::numTruncated() const
{
    return static_cast<index_t>(std::count(m_truncated.begin(), m_truncated.end(), true));
}

template<class T>
index_t gsMPBES<T>::numFunctionsWithComponents(index_t n) const
{
    index_t c = 0;
    for (size_t g = 0; g != m_components.size(); ++g)
        if ( (n >= 3 && static_cast<index_t>(m_components[g].size()) >= 3) ||
             static_cast<index_t>(m_components[g].size()) == n ) ++c;
    return c;
}

template<class T>
gsMatrix<T> gsMPBES<T>::representation(const gsMultiPatch<T> & mp, T * residual) const
{
    const size_t nP = m_thb.nBases();
    GISMO_ENSURE(mp.nPatches() == nP, "gsMPBES: the multipatch has "<<mp.nPatches()
                 <<" patches, the basis "<<nP<<".");
    index_t nLocal = 0;
    for (size_t p = 0; p != nP; ++p) nLocal += m_levels[p][maxLevel(p)].size();
    gsMatrix<T> local(nLocal, mp.geoDim());
    index_t off = 0;
    for (size_t p = 0; p != nP; ++p)
    {
        const index_t n = m_levels[p][maxLevel(p)].size();
        gsMatrix<T> coefs;
        if ( const gsTHBSpline<2,T> * thb = dynamic_cast<const gsTHBSpline<2,T>*>(&mp.patch(p)) )
        {
            gsTHBSpline<2,T> copy(*thb);
            gsTensorBSpline<2,T> fine;
            copy.convertToBSpline(fine);
            coefs = fine.coefs();
        }
        else
            coefs = mp.patch(p).coefs();
        GISMO_ENSURE(coefs.rows() == n, "gsMPBES: patch "<<p<<" of the multipatch has "<<coefs.rows()
                     <<" coefficients on its finest level, the basis expects "<<n<<".");
        local.middleRows(off, n) = coefs;
        off += n;
    }
    gsMatrix<T> global;
    Base::local_coef_to_global_coef(local, global);
    if (residual)
    {
        gsMatrix<T> back;
        Base::global_coef_to_local_coef(global, back);
        const T nrm = local.norm();
        *residual = (back - local).norm() / (nrm > 0 ? nrm : (T)1);
    }
    return global;
}

template<class T>
T gsMPBES<T>::partitionOfUnityDefect(index_t samples) const
{
    T defect = 0;
    gsMatrix<T> vals;
    for (size_t p = 0; p != m_thb.nBases(); ++p)
    {
        const gsMatrix<T> supp = m_levels[p][0].support();
        gsVector<unsigned> np(2);
        np.setConstant(samples);
        const gsMatrix<T> pts = gsPointGrid<T>(supp.col(0), supp.col(1), np);
        Base::eval_into(static_cast<index_t>(p), pts, vals);
        for (index_t i = 0; i != vals.cols(); ++i)
            defect = math::max(defect, math::abs(vals.col(i).sum() - (T)1));
    }
    return defect;
}

template<class T>
T gsMPBES<T>::interfaceJump(index_t samples, int seed) const
{
    std::mt19937 gen(seed);
    std::uniform_real_distribution<T> dist(-1, 1);
    gsMatrix<T> coefs(this->size(), 1);
    for (index_t i = 0; i != coefs.rows(); ++i) coefs(i,0) = dist(gen);
    gsMappedSpline<2,T> f(*this, coefs);

    T jump = 0;
    const gsBoxTopology & topol = m_thb.topology();
    gsMatrix<T> u1, u2, v1, v2;
    for (gsBoxTopology::const_iiterator it = topol.iBegin(); it != topol.iEnd(); ++it)
    {
        const boundaryInterface & bi = *it;
        const index_t p1 = bi.first().patch, p2 = bi.second().patch;
        const gsMatrix<T> box1 = m_levels[p1][0].support(), box2 = m_levels[p2][0].support();
        // points on the side of patch 1
        const short_t dir = bi.first().side().direction();
        const T fixed = bi.first().side().parameter() ? box1(dir,1) : box1(dir,0);
        u1.resize(2, samples);
        for (index_t i = 0; i != samples; ++i)
        {
            const T t = (T)(i + 1) / (T)(samples + 1);
            u1(dir,   i) = fixed;
            u1(1-dir, i) = box1(1-dir,0) + t * (box1(1-dir,1) - box1(1-dir,0));
        }
        gsAffineFunction<T> amap(bi.dirMap(), bi.dirOrientation(), box1, box2);
        amap.eval_into(u1, u2);
        f.piece(p1).eval_into(u1, v1);
        f.piece(p2).eval_into(u2, v2);
        jump = math::max(jump, (v1 - v2).cwiseAbs().maxCoeff());
    }
    return jump;
}

template<class T>
std::ostream & gsMPBES<T>::print(std::ostream & os) const
{
    os << "gsMPBES: " << this->size() << " functions on " << m_thb.nBases() << " patches, "
       << m_numLevels << " levels (per level:";
    for (index_t l = 0; l != m_numLevels; ++l) os << " " << numFunctionsOfLevel(l);
    os << "), " << numFunctionsWithComponents(1) << " interior, " << numFunctionsWithComponents(2)
       << " interface, " << numFunctionsWithComponents(3) << " vertex functions, "
       << numTruncated() << " truncated\n";
    return os;
}

} // namespace gismo
