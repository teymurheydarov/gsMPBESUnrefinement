/** @file gsMPBESUnrefinement.hpp

    @brief Implementation of gsMPBESUnrefinement

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): T. Heydarov
*/

#pragma once

#include <gsMPBESUnrefinement/src/gsMPBESUnrefinement.h>
#include <gsMPBESUnrefinement/src/gsMPBESSampler.h>
#include <gsHSplines/gsTHBSplineBasis.h>
#include <gsHSplines/gsTHBSpline.h>
#include <gsDomain/gsDomain.h>
#include <gsAssembler/gsGaussRule.h>
#include <gsAssembler/gsLobattoRule.h>
#include <gsMatrix/gsSparseSolver.h>
#include <gsUtils/gsStopwatch.h>
#include <unordered_map>
#include <set>
#include <cmath>
#include <sstream>
#include <iomanip>

namespace gismo
{

namespace internal
{
    inline long long componentKey(index_t patch, index_t level, index_t index)
    {
        return (static_cast<long long>(patch) << 40) | (static_cast<long long>(level) << 32)
             | static_cast<long long>(static_cast<unsigned int>(index));
    }

    /// positive-area overlap of two boxes [u1,u2]x[v1,v2]
    template<class T>
    bool boxesOverlap(const std::array<T,4> & a, const std::array<T,4> & b)
    {
        return math::min(a[1], b[1]) > math::max(a[0], b[0]) &&
               math::min(a[3], b[3]) > math::max(a[2], b[2]);
    }

    template<class T>
    std::array<T,4> supportBox(const gsTensorBSplineBasis<2,T> & b, index_t i)
    {
        const gsMatrix<T> s = b.support(i);
        return { s(0,0), s(0,1), s(1,0), s(1,1) };
    }
}

template<class T>
gsMPBESUnrefinement<T>::gsMPBESUnrefinement(const gsMultiPatch<T> & p0, const Options & opt)
: m_p0(p0), m_opt(opt)
{
    GISMO_ENSURE(m_p0.nPatches() > 0, "gsMPBESUnrefinement: no patches.");
    if (m_p0.nInterfaces() == 0 && m_p0.nPatches() > 1)
        m_p0.computeTopology();
    if (m_p0.nBoundary() == 0 && m_p0.nInterfaces() == 0)
        m_p0.computeTopology();
    const size_t nP = m_p0.nPatches();

    // the feature set
    if (m_opt.feature.empty())
        m_feature = m_p0.boundaries();
    else
    {
        for (const patchSide & ps : m_opt.feature)
        {
            GISMO_ENSURE(ps.patch >= 0 && ps.patch < static_cast<index_t>(nP),
                         "gsMPBESUnrefinement: feature side on patch " << ps.patch << " of " << nP << ".");
            bool interface = false;
            for (const boundaryInterface & bi : m_p0.interfaces())
                if (bi.first() == ps || bi.second() == ps) { interface = true; break; }
            if (interface)
            {
                gsWarn << "gsMPBESUnrefinement: the feature side " << ps << " is an interface and is skipped.\n";
                continue;
            }
            bool dup = false;
            for (const patchSide & q : m_feature) if (q == ps) dup = true;
            if (!dup) m_feature.push_back(ps);
        }
    }
    if (m_opt.interfacesInFeature)
        for (const boundaryInterface & bi : m_p0.interfaces())
            m_feature.push_back(bi.first());

    // the grids of (C2) and of the fit, and the points of (C3)
    m_gridSize.resize(nP);
    m_grid.resize(nP);
    m_gridTarget.resize(nP);
    m_featurePts.assign(nP, gsMatrix<T>(2, 0));
    m_featureTarget.assign(nP, gsMatrix<T>(0, 2));
    for (size_t p = 0; p != nP; ++p)
    {
        const T n = static_cast<T>(m_p0.patch(p).basis().size());
        m_gridSize[p] = math::max(index_t(5), static_cast<index_t>(math::ceil(m_opt.sampleFactor * math::sqrt(n))));
        m_grid[p] = gsMPBESSampler<T>::uniformGrid(m_gridSize[p], m_gridSize[p]);
        gsMatrix<T> xy;
        m_p0.patch(p).eval_into(m_grid[p], xy);
        m_gridTarget[p] = xy.transpose();
    }
    for (const patchSide & ps : m_feature)
    {
        const gsMatrix<T> pts = gsMPBESSampler<T>::sidePoints(ps.side(), m_opt.featurePoints);
        gsMatrix<T> & all = m_featurePts[ps.patch];
        const index_t c = all.cols();
        all.conservativeResize(2, c + pts.cols());
        all.rightCols(pts.cols()) = pts;
        gsMatrix<T> xy;
        m_p0.patch(ps.patch).eval_into(pts, xy);
        gsMatrix<T> & tg = m_featureTarget[ps.patch];
        tg.conservativeResize(c + pts.cols(), 2);
        tg.bottomRows(pts.cols()) = xy.transpose();
    }
}

template<class T>
void gsMPBESUnrefinement<T>::log(int level, const std::string & line) const
{
    if (m_opt.verbosity >= level)
        gsInfo << line << "\n";
    if (m_opt.log)
        (*m_opt.log) << line << "\n";
}

template<class T>
bool gsMPBESUnrefinement<T>::lexicographicOrder(const Cell & a, const Cell & b)
{
    if (a.patch != b.patch) return a.patch < b.patch;
    if (a.y1 != b.y1) return a.y1 > b.y1;
    return a.x1 > b.x1;
}

template<class T>
std::array<T,4> gsMPBESUnrefinement<T>::cellBox(const gsMPBES<T> & basis, const Cell & c)
{
    const gsTensorBSplineBasis<2,T> & lev = basis.tensorLevel(c.patch, c.level);
    const T nu = static_cast<T>(lev.knots(0).numElements()), nv = static_cast<T>(lev.knots(1).numElements());
    return { T(c.x1) / nu, T(c.x1 + 1) / nu, T(c.y1) / nv, T(c.y1 + 1) / nv };
}

template<class T>
std::vector<typename gsMPBESUnrefinement<T>::Cell>
gsMPBESUnrefinement<T>::eligibleCells(const gsMPBES<T> & basis, index_t level)
{
    std::vector<Cell> cells;
    for (size_t p = 0; p != basis.nPatches(); ++p)
    {
        if (basis.maxLevel(static_cast<index_t>(p)) < level + 1) continue;
        const gsTensorBSplineBasis<2,T> & lev = basis.tensorLevel(static_cast<index_t>(p), level);
        const index_t nu = static_cast<index_t>(lev.knots(0).numElements());
        const index_t nv = static_cast<index_t>(lev.knots(1).numElements());
        const gsHTree<2,index_t> & tree = basis.thbBasis(static_cast<index_t>(p)).tree();
        gsVector<index_t,2> lo, up;
        for (index_t y = nv - 1; y >= 0; --y)
            for (index_t x = nu - 1; x >= 0; --x)
            {
                lo << x, y;
                up << x + 1, y + 1;
                // c in Omega^{level+1} and the interior of c not meeting Omega^{level+2}:
                // the lowest and the highest level inside c are both level + 1
                if (tree.query3(lo, up, static_cast<int>(level)) == static_cast<int>(level) + 1 &&
                    tree.query4(lo, up, static_cast<int>(level)) == static_cast<int>(level) + 1)
                    cells.push_back({static_cast<index_t>(p), level, x, y});
            }
    }
    return cells;
}

template<class T>
gsMultiBasis<T> gsMPBESUnrefinement<T>::coarsen(const gsMultiBasis<T> & thb, const std::vector<Cell> & group)
{
    gsMultiBasis<T> result(thb);
    std::map< index_t, std::vector<index_t> > boxes;    // patch -> boxes (level, x1, y1, x2, y2)
    for (const Cell & c : group)
    {
        std::vector<index_t> & b = boxes[c.patch];
        b.push_back(c.level);
        b.push_back(c.x1); b.push_back(c.y1);
        b.push_back(c.x1 + 1); b.push_back(c.y1 + 1);
    }
    for (typename std::map< index_t, std::vector<index_t> >::const_iterator it = boxes.begin(); it != boxes.end(); ++it)
    {
        gsHTensorBasis<2,T> * h = dynamic_cast<gsHTensorBasis<2,T>*>(&result.basis(it->first));
        GISMO_ENSURE(h, "gsMPBESUnrefinement: patch " << it->first << " is not hierarchical.");
        h->unrefineElements(it->second);
    }
    return result;
}

template<class T>
gsMatrix<T> gsMPBESUnrefinement<T>::transferCoefficients(const gsMPBES<T> & from, const gsMatrix<T> & coefs,
                                                         const gsMPBES<T> & to, std::vector<index_t> * newFunctions)
{
    GISMO_ASSERT(coefs.rows() == from.size(), "gsMPBESUnrefinement: coefficients of the source basis expected.");
    std::unordered_map<long long, index_t> lookup;
    for (index_t g = 0; g != from.size(); ++g)
        for (const typename gsMPBES<T>::Component & c : from.components(g))
            lookup[internal::componentKey(c.patch, c.level, c.index)] = g;
    gsMatrix<T> result;
    result.setZero(to.size(), coefs.cols());
    if (newFunctions) newFunctions->clear();
    for (index_t g = 0; g != to.size(); ++g)
    {
        bool found = false;
        for (const typename gsMPBES<T>::Component & c : to.components(g))
        {
            std::unordered_map<long long, index_t>::const_iterator it = lookup.find(internal::componentKey(c.patch, c.level, c.index));
            if (it != lookup.end())
            {
                result.row(g) = coefs.row(it->second);
                found = true;
                break;
            }
        }
        if (!found && newFunctions)
            newFunctions->push_back(g);
    }
    return result;
}

template<class T>
gsMultiPatch<T> gsMPBESUnrefinement<T>::thbPatches(const gsMPBES<T> & basis, const gsMatrix<T> & coefs)
{
    // every patch of the parameterisation lies in the THB-spline space
    // of its hierarchy: least-squares collocation at (d+1)^2 points of
    // every element, which is exact
    gsMultiPatch<T> result;
    const std::vector< gsMatrix<T> > pts = elementPoints(basis, std::vector< std::vector< std::array<T,4> > >());
    for (size_t p = 0; p != basis.nPatches(); ++p)
    {
        const gsTHBSplineBasis<2,T> & thb = basis.thbBasis(static_cast<index_t>(p));
        const gsMatrix<T> & u = pts[p];
        const index_t n = u.cols();
        // the THB basis at the points
        gsMatrix<index_t> act;
        gsMatrix<T> val;
        thb.active_into(u, act);
        thb.eval_into(u, val);
        std::vector< gsEigen::Triplet<T> > trip;
        for (index_t k = 0; k != n; ++k)
            for (index_t i = 0; i != act.rows(); ++i)
                if (val(i, k) != 0)
                    trip.push_back(gsEigen::Triplet<T>(k, act(i, k), val(i, k)));
        gsSparseMatrix<T> A(n, thb.size());
        A.setFromTriplets(trip.begin(), trip.end());
        const gsMatrix<T> rhs = gsMPBESSampler<T>::evalGeometry(basis, static_cast<index_t>(p), u, coefs).transpose();
        gsSparseMatrix<T> AtA = A.transpose() * A;
        typename gsSparseSolver<T>::SimplicialLDLT solver(AtA);
        gsMatrix<T> x = solver.solve(A.transpose() * rhs);
        GISMO_ENSURE(solver.info() == gsEigen::Success && x.allFinite(),
                     "gsMPBESUnrefinement: the THB-spline representation of patch " << p << " failed.");
        result.addPatch(gsTHBSpline<2,T>(thb, x));
    }
    result.computeTopology();
    return result;
}

template<class T>
typename gsMPBESUnrefinement<T>::Check
gsMPBESUnrefinement<T>::check(const gsMPBES<T> & basis, const gsMatrix<T> & coefs, bool certify) const
{
    Check ch;
    const size_t nP = m_p0.nPatches();
    for (size_t p = 0; p != nP; ++p)
    {
        const gsMatrix<T> xy = gsMPBESSampler<T>::evalGeometry(basis, static_cast<index_t>(p), m_grid[p], coefs).transpose();
        ch.globalError = math::max(ch.globalError, (xy - m_gridTarget[p]).rowwise().norm().maxCoeff());
        if (m_featurePts[p].cols() > 0)
        {
            const gsMatrix<T> f = gsMPBESSampler<T>::evalGeometry(basis, static_cast<index_t>(p), m_featurePts[p], coefs).transpose();
            ch.featureError = math::max(ch.featureError, (f - m_featureTarget[p]).rowwise().norm().maxCoeff());
        }
    }
    if (!(ch.globalError == ch.globalError) || !(ch.featureError == ch.featureError))   // NaN
    {
        ch.globalError = ch.featureError = std::numeric_limits<T>::infinity();
        return ch;
    }
    const bool accurate = ch.globalError <= m_opt.epsilonGlobal && ch.featureError <= m_opt.epsilonFeature;
    if (certify)
    {
        gsBernsteinCertificate<T> cert(m_opt.certificateDepth);
        typename gsBernsteinCertificate<T>::Result r = cert.certify(basis, coefs, basis.thbBases());
        ch.certified = r.certified;
        ch.notCertified = r.notCertified;
    }
    ch.valid = accurate && ch.certified;
    return ch;
}

template<class T>
std::vector<typename gsMPBESUnrefinement<T>::Cell>
gsMPBESUnrefinement<T>::grendasSet(const gsMPBES<T> & basis, const gsMatrix<T> & coefs,
                                   const std::vector<Cell> & working) const
{
    if (working.empty()) return std::vector<Cell>();
    // eta_c = int_c |J|_F^2 with a 3 x 3 Gauss rule, evaluated patch by patch
    gsVector<index_t> numNodes(2);
    numNodes.setConstant(3);
    gsGaussRule<T> rule(numNodes);
    std::vector<T> eta(working.size(), T(0));
    std::map< index_t, std::vector<size_t> > byPatch;
    for (size_t k = 0; k != working.size(); ++k) byPatch[working[k].patch].push_back(k);
    gsMatrix<T> nodes;
    gsVector<T> w;
    gsVector<T> lo(2), up(2);
    for (typename std::map< index_t, std::vector<size_t> >::const_iterator it = byPatch.begin(); it != byPatch.end(); ++it)
    {
        const std::vector<size_t> & ks = it->second;
        gsMatrix<T> all(2, 9 * ks.size());
        std::vector<T> weights(9 * ks.size());
        for (size_t j = 0; j != ks.size(); ++j)
        {
            const std::array<T,4> b = cellBox(basis, working[ks[j]]);
            lo << b[0], b[2];
            up << b[1], b[3];
            rule.mapTo(lo, up, nodes, w);
            all.middleCols(9 * j, 9) = nodes;
            for (index_t q = 0; q != 9; ++q) weights[9 * j + q] = w[q];
        }
        std::vector< gsSparseMatrix<T> > B;
        gsMPBESSampler<T>::sample(basis, it->first, all, 1, B);
        const gsMatrix<T> Ju = B[1] * coefs, Jv = B[2] * coefs;
        for (size_t j = 0; j != ks.size(); ++j)
            for (index_t q = 0; q != 9; ++q)
                eta[ks[j]] += weights[9 * j + q] * (Ju.row(9 * j + q).squaredNorm() + Jv.row(9 * j + q).squaredNorm());
    }
    T total = 0;
    for (size_t k = 0; k != eta.size(); ++k) total += eta[k];
    std::vector<size_t> order(working.size());
    for (size_t k = 0; k != order.size(); ++k) order[k] = k;
    std::stable_sort(order.begin(), order.end(), [&eta](size_t a, size_t b) { return eta[a] < eta[b]; });
    std::vector<Cell> group;
    T running = 0;
    for (size_t k = 0; k != order.size(); ++k)
    {
        group.push_back(working[order[k]]);
        running += eta[order[k]];
        if (running >= m_opt.grendasDelta * total) break;
    }
    return group;
}

template<class T>
std::vector< gsMatrix<T> > gsMPBESUnrefinement<T>::elementPoints(const gsMPBES<T> & basis,
                                                               const std::vector< std::vector< std::array<T,4> > > & domain)
{
    const size_t nP = basis.nPatches();
    std::vector< gsMatrix<T> > uv(nP, gsMatrix<T>(2, 0));
    for (size_t p = 0; p != nP; ++p)
    {
        if (!domain.empty() && domain[p].empty()) continue;
        const gsTHBSplineBasis<2,T> & thb = basis.thbBasis(static_cast<index_t>(p));
        const index_t du = thb.degree(0) + 1, dv = thb.degree(1) + 1;
        std::vector< gsMatrix<T> > pts;
        index_t n = 0;
        typename gsDomain<T>::Ptr dom = thb.domain();
        for (auto it = dom->beginAll(); it != dom->endAll(); ++it)
        {
            const gsVector<T> lo = it.lowerCorner(), up = it.upperCorner();
            if (!domain.empty())
            {
                const std::array<T,4> e = { lo(0), up(0), lo(1), up(1) };
                bool meets = false;
                for (const std::array<T,4> & b : domain[p])
                    if (internal::boxesOverlap(e, b)) { meets = true; break; }
                if (!meets) continue;
            }
            gsMatrix<T> g(2, du * dv);
            for (index_t j = 0; j != dv; ++j)
                for (index_t i = 0; i != du; ++i)
                {
                    g(0, j * du + i) = lo(0) + (up(0) - lo(0)) * T(2 * i + 1) / T(2 * du);
                    g(1, j * du + i) = lo(1) + (up(1) - lo(1)) * T(2 * j + 1) / T(2 * dv);
                }
            pts.push_back(g);
            n += g.cols();
        }
        uv[p].resize(2, n);
        index_t c = 0;
        for (const gsMatrix<T> & g : pts) { uv[p].middleCols(c, g.cols()) = g; c += g.cols(); }
    }
    return uv;
}

template<class T>
void gsMPBESUnrefinement<T>::gaussSamples(const gsMPBES<T> & basis,
                                          const std::vector< std::vector< std::array<T,4> > > & domain,
                                          std::vector< gsMatrix<T> > & uv, std::vector< gsVector<T> > & w) const
{
    const size_t nP = basis.nPatches();
    uv.assign(nP, gsMatrix<T>(2, 0));
    w.assign(nP, gsVector<T>());
    for (size_t p = 0; p != nP; ++p)
    {
        if (!domain.empty() && domain[p].empty()) continue;
        const gsTHBSplineBasis<2,T> & thb = basis.thbBasis(static_cast<index_t>(p));
        gsVector<index_t> numNodes(2);
        numNodes[0] = (m_opt.gaussNodes > 0) ? m_opt.gaussNodes : thb.degree(0) + 1;
        numNodes[1] = (m_opt.gaussNodes > 0) ? m_opt.gaussNodes : thb.degree(1) + 1;
        // Gauss-Lobatto: the points on an edge shared by two elements enter twice, each
        // with the weight of its element, which is the composite rule
        typename gsQuadRule<T>::uPtr rule;
        if (m_opt.samples == GaussLobattoPoints) rule = gsLobattoRule<T>::make(numNodes);
        else                                     rule = gsGaussRule<T>::make(numNodes);
        std::vector< gsMatrix<T> > pts;
        std::vector<T> weights;
        gsMatrix<T> nodes;
        gsVector<T> ew;
        typename gsDomain<T>::Ptr dom = thb.domain();
        for (auto it = dom->beginAll(); it != dom->endAll(); ++it)
        {
            const gsVector<T> lo = it.lowerCorner(), up = it.upperCorner();
            if (!domain.empty())
            {
                const std::array<T,4> e = { lo(0), up(0), lo(1), up(1) };
                bool meets = false;
                for (const std::array<T,4> & b : domain[p])
                    if (internal::boxesOverlap(e, b)) { meets = true; break; }
                if (!meets) continue;
            }
            rule->mapTo(lo, up, nodes, ew);
            pts.push_back(nodes);
            for (index_t k = 0; k != ew.size(); ++k) weights.push_back(ew[k]);
        }
        index_t n = 0;
        for (const gsMatrix<T> & g : pts) n += g.cols();
        uv[p].resize(2, n);
        index_t c = 0;
        for (const gsMatrix<T> & g : pts) { uv[p].middleCols(c, g.cols()) = g; c += g.cols(); }
        if (!m_opt.unitWeights)
        {
            w[p] = gsAsConstVector<T>(weights);
            continue;
        }
        // every distinct point once with weight 1: the Gauss-Lobatto points on an edge
        // shared by two elements (and at a corner of four) are merged
        std::set< std::pair<long long, long long> > seen;
        index_t m = 0;
        for (index_t k = 0; k != n; ++k)
        {
            const std::pair<long long, long long> key(std::llround(uv[p](0, k) * T(1e10)),
                                                      std::llround(uv[p](1, k) * T(1e10)));
            if (seen.insert(key).second) uv[p].col(m++) = uv[p].col(k);
        }
        uv[p].conservativeResize(2, m);
        w[p].setOnes(m);
    }
}

template<class T>
void gsMPBESUnrefinement<T>::fittingSetup(const gsMPBES<T> & candidate, const std::vector<Cell> & group,
                                          const std::vector<index_t> & newFunctions,
                                          gsMPBESProjection<T> & proj, index_t & numFree) const
{
    const size_t nP = m_p0.nPatches();
    const index_t N = candidate.size();
    const std::vector< std::vector< std::array<T,4> > > whole;
    const bool gauss = (m_opt.samples != UniformGrid);
    auto setGlobalSamples = [&]()
    {
        if (gauss)
        {
            std::vector< gsMatrix<T> > uv;
            std::vector< gsVector<T> > w;
            gaussSamples(candidate, whole, uv, w);
            proj.setSamples(uv, w);
        }
        else
            proj.setSamples(m_grid);
    };
    if (m_opt.lambda < 0)
    {
        proj.setFreeFunctions(std::vector<bool>());
        setGlobalSamples();
        proj.setAugmentation([&candidate, whole]() { return elementPoints(candidate, whole); });
        proj.setFittingDomain(whole);
        numFree = N;
        return;
    }
    // the extended cells E(c)
    std::vector< std::vector< std::array<T,4> > > ext(nP);
    for (const Cell & c : group)
    {
        const std::array<T,4> b = cellBox(candidate, c);
        const T du = b[1] - b[0], dv = b[3] - b[2], l = static_cast<T>(m_opt.lambda);
        ext[c.patch].push_back({ math::max(T(0), b[0] - l * du), math::min(T(1), b[1] + l * du),
                                 math::max(T(0), b[2] - l * dv), math::min(T(1), b[3] + l * dv) });
    }
    // the interior set: functions with a component on the cell's patch meeting E(c), and the new functions
    std::vector<bool> free(N, false);
    for (index_t g = 0; g != N; ++g)
        for (const typename gsMPBES<T>::Component & c : candidate.components(g))
        {
            if (ext[c.patch].empty()) continue;
            const std::array<T,4> s = internal::supportBox(candidate.tensorLevel(c.patch, c.level), c.index);
            bool meets = false;
            for (const std::array<T,4> & e : ext[c.patch])
                if (internal::boxesOverlap(s, e)) { meets = true; break; }
            if (meets) { free[g] = true; break; }
        }
    for (index_t g : newFunctions) free[g] = true;
    numFree = 0;
    for (index_t g = 0; g != N; ++g) if (free[g]) ++numFree;
    if (numFree == N)
    {
        proj.setFreeFunctions(std::vector<bool>());
        setGlobalSamples();
        proj.setAugmentation([&candidate, whole]() { return elementPoints(candidate, whole); });
        proj.setFittingDomain(whole);
        return;
    }
    // the fitting domain: the supports of the interior functions on every patch
    std::vector< std::vector< std::array<T,4> > > domain(nP);
    for (index_t g = 0; g != N; ++g)
        if (free[g])
            for (const typename gsMPBES<T>::Component & c : candidate.components(g))
                domain[c.patch].push_back(internal::supportBox(candidate.tensorLevel(c.patch, c.level), c.index));
    // the sample points: on every patch of the domain, the uniform grid of the box covering
    // Delta(c)|p at the density of the grid of (C2) (q_p points per unit length, at least 3)
    std::vector< gsMatrix<T> > uv(nP, gsMatrix<T>(2, 0));
    for (size_t p = 0; p != nP; ++p)
    {
        if (domain[p].empty()) continue;
        T u0 = 1, u1 = 0, v0 = 1, v1 = 0;
        for (const std::array<T,4> & b : domain[p])
        {
            u0 = math::min(u0, b[0]); u1 = math::max(u1, b[1]);
            v0 = math::min(v0, b[2]); v1 = math::max(v1, b[3]);
        }
        const index_t k = math::max(index_t(3), m_gridSize[p]);
        const index_t ku = math::max(index_t(3), static_cast<index_t>(math::ceil(k * (u1 - u0) - 1e-12)));
        const index_t kv = math::max(index_t(3), static_cast<index_t>(math::ceil(k * (v1 - v0) - 1e-12)));
        uv[p] = gsMPBESSampler<T>::uniformGrid(ku, kv, u0, u1, v0, v1);
    }
    proj.setFreeFunctions(free);
    if (gauss)
    {
        std::vector< gsMatrix<T> > guv;
        std::vector< gsVector<T> > gw;
        gaussSamples(candidate, domain, guv, gw);
        proj.setSamples(guv, gw);
    }
    else
        proj.setSamples(uv);
    proj.setAugmentation([&candidate, domain]() { return elementPoints(candidate, domain); });
    proj.setFittingDomain(domain);
}

template<class T>
bool gsMPBESUnrefinement<T>::project(const gsMPBES<T> & candidate, gsMatrix<T> & coefs,
                                     const std::vector<Cell> & group, const std::vector<index_t> & newFunctions,
                                     Record & rec) const
{
    gsMPBESProjection<T> proj(candidate, m_p0, m_opt.projection);
    index_t numFree = 0;
    gsStopwatch clock;
    fittingSetup(candidate, group, newFunctions, proj, numFree);
    rec.secondsAssembly += clock.stop();
    rec.numFree = numFree;
    const size_t nP = m_p0.nPatches();
    auto timedCheck = [&](const gsMatrix<T> & c) -> Check
    {
        gsStopwatch w;
        const Check ch = check(candidate, c, true);
        rec.secondsCheck += w.stop();
        return ch;
    };
    auto rho = [&](const Check & ch) -> T
    {
        // rho = 100 n_- / n_J with n_J = nP m^2 the points of a uniform m x m regularity grid per patch
        const index_t m = math::max(index_t(20), static_cast<index_t>(math::ceil(2 * math::sqrt(T(candidate.size()) / T(nP)))));
        return T(100) * T(ch.notCertified) / T(nP * m * m);
    };
    std::ostringstream detail;

    // Projection[1]: least squares
    clock.restart();
    proj.leastSquares(coefs);
    rec.secondsSolve += clock.stop();
    rec.numSamples = proj.numSamples();
    rec.augmented = proj.augmentedSamples();
    rec.stage = 1;
    rec.check = timedCheck(coefs);
    detail << "LS: e_g " << rec.check.globalError << ", e_f " << rec.check.featureError << ", n_- " << rec.check.notCertified;
    if (rec.check.valid) { log(2, "    " + detail.str()); return true; }
    if (rec.check.certified || m_opt.stages < 2)
    {
        // regular but outside the tolerances: the fit is optimal for the accuracy, the
        // optimisation stages cannot help
        log(2, "    " + detail.str());
        return false;
    }
    const gsMatrix<T> lsCoefs = coefs;

    // Projection[2]: linear optimisation
    typename gsMPBESProjection<T>::Weights w = gsMPBESProjection<T>::adaptiveWeights(
        rec.check.globalError, rec.check.featureError, rho(rec.check), rec.check.notCertified,
        m_opt.projection.bLinear, true);
    clock.restart();
    proj.linearOptimisation(coefs, w);
    rec.secondsOptimisation += clock.stop();
    rec.stage = 2;
    rec.check = timedCheck(coefs);
    detail << " | LO (w_u " << w.uniformity << ", w_l " << w.length << "): e_g " << rec.check.globalError
           << ", e_f " << rec.check.featureError << ", n_- " << rec.check.notCertified;
    if (rec.check.valid) { log(2, "    " + detail.str()); return true; }
    if (rec.check.certified || m_opt.stages < 3)
    {
        log(2, "    " + detail.str());
        return false;
    }

    // Projection[3]: non-linear optimisation, from the least-squares candidate with the weights
    // determined by the linear-stage candidate
    w = gsMPBESProjection<T>::adaptiveWeights(
        rec.check.globalError, rec.check.featureError, rho(rec.check), rec.check.notCertified,
        m_opt.projection.bNonlinear, false);
    coefs = lsCoefs;
    Check last;
    typename gsMPBESProjection<T>::Validator validator = [&](const gsMatrix<T> & c) -> bool
    {
        last = timedCheck(c);
        return last.valid;
    };
    clock.restart();
    const double checksBefore = rec.secondsCheck;
    typename gsMPBESProjection<T>::Stats st = proj.nonlinearOptimisation(coefs, w, validator);
    rec.secondsOptimisation += clock.stop() - (rec.secondsCheck - checksBefore);
    rec.stage = 3;
    rec.check = (st.iterations > 0) ? last : timedCheck(coefs);
    detail << " | NLO (w_u " << w.uniformity << ", w_l " << w.length << ", w_o " << w.orthogonality
           << ", w_s " << w.skewness << ", w_e " << w.eccentricity << "; " << st.iterations << " it., "
           << st.stop << "): e_g " << rec.check.globalError << ", e_f " << rec.check.featureError
           << ", n_- " << rec.check.notCertified;
    log(2, "    " + detail.str());
    return rec.check.valid;
}

template<class T>
typename gsMPBESUnrefinement<T>::Result gsMPBESUnrefinement<T>::run()
{
    gsStopwatch clock;
    Result res;
    gsMPBES<T>  basis(m_p0);
    gsMatrix<T> coefs = basis.representation(m_p0);
    res.initialSize = basis.size();
    {
        const Check c0 = check(basis, coefs, true);
        std::ostringstream s;
        s << "gsMPBESUnrefinement: " << m_p0.nPatches() << " patches, " << basis.numLevels() << " levels, "
          << basis.size() << " MPBES functions; p0 in the MPBES space: e_g " << c0.globalError
          << ", e_f " << c0.featureError << ", " << (c0.certified ? "certified" : "NOT certified")
          << " (n_- " << c0.notCertified << ")";
        log(1, s.str());
        s.str("");
        s << "  epsilon_g " << m_opt.epsilonGlobal << ", epsilon_f " << m_opt.epsilonFeature
          << ", feature " << m_feature.size() << " sides, selection " << (m_opt.selection == Grendas ? "Grendas" : "lexicographic")
          << ", lambda " << m_opt.lambda << ", stages " << m_opt.stages << ", solver " << m_opt.projection.solver
          << ", fit samples " << (m_opt.samples == GaussPoints ? "Gauss points"
                                  : m_opt.samples == GaussLobattoPoints ? "Gauss-Lobatto points" : "uniform grid")
          << (m_opt.samples != UniformGrid && m_opt.unitWeights ? " (unit weights)" : "")
          << ", functional quadrature " << m_opt.projection.quadratureNodes << " nodes per direction";
        log(1, s.str());
    }

    const index_t numLevels = basis.numLevels();
    for (index_t level = numLevels - 2; level >= 0; --level)
    {
        // Algorithm 1, lines 10-28: the eligible set E of the level, passes over it while
        // the preceding pass accepted a coarsening; the working set of a pass starts as E,
        // a group leaves it when selected, a rejected group of several cells returns to it,
        // and E becomes the working set at every acceptance
        std::vector<Cell> eligible = eligibleCells(basis, level);
        {
            std::ostringstream s;
            s << "level " << level + 1 << " (cells of level " << level << " refined to " << level + 1 << "): "
              << eligible.size() << " eligible cells";
            log(1, s.str());
        }
        bool accepted = true;
        index_t pass = 0;
        while (!eligible.empty() && accepted)
        {
            ++pass;
            if (pass > 1)
            {
                std::ostringstream s;
                s << "  pass " << pass << " over the " << eligible.size() << " remaining eligible cells";
                log(1, s.str());
            }
            std::vector<Cell> working = eligible;
            bool grendas = (m_opt.selection == Grendas);
            accepted = false;
            while (!working.empty())
            {
                std::vector<Cell> group;
                bool viaGrendas = false;
                if (grendas)
                {
                    gsStopwatch w;
                    group = grendasSet(basis, coefs, working);
                    res.timing.selection += w.stop();
                    viaGrendas = !group.empty();
                }
                if (group.empty())
                {
                    std::sort(working.begin(), working.end(), lexicographicOrder);
                    group.push_back(working.front());
                }
                grendas = false;
                // G leaves the working set
                {
                    std::vector<Cell> rest;
                    for (const Cell & c : working)
                    {
                        bool inGroup = false;
                        for (const Cell & g : group)
                            if (g.patch == c.patch && g.x1 == c.x1 && g.y1 == c.y1) { inGroup = true; break; }
                        if (!inGroup) rest.push_back(c);
                    }
                    working.swap(rest);
                }
                gsStopwatch candidateClock;
                Record rec;
                rec.level = level;
                rec.pass = pass;
                rec.group = group;
                rec.viaGrendas = viaGrendas;
                // the candidate hierarchy and basis
                gsMPBES<T> candidate(coarsen(basis.thbBases(), group));
                rec.size = candidate.size();
                std::vector<index_t> newFunctions;
                gsMatrix<T> candCoefs = transferCoefficients(basis, coefs, candidate, &newFunctions);
                rec.secondsBasis = candidateClock.stop();
                const bool valid = project(candidate, candCoefs, group, newFunctions, rec);
                rec.accepted = valid;
                rec.seconds = candidateClock.stop();
                res.timing.basis += rec.secondsBasis;
                res.timing.assembly += rec.secondsAssembly;
                res.timing.solve += rec.secondsSolve;
                res.timing.check += rec.secondsCheck;
                res.timing.optimisation += rec.secondsOptimisation;
                {
                    std::set<index_t> patches;
                    for (const Cell & c : group) patches.insert(c.patch);
                    std::ostringstream s;
                    s << "  " << (viaGrendas ? "Grendas" : "lex") << " group of " << group.size() << " cell(s) on "
                      << patches.size() << " patch(es) (first: patch " << group.front().patch << ", cell ("
                      << group.front().x1 << "," << group.front().y1 << ")): " << candidate.size() << " functions, "
                      << rec.numFree << " fitted, " << rec.numSamples << (rec.augmented ? " samples (grid + element points); stage " : " samples; stage ") << rec.stage << ": e_g "
                      << rec.check.globalError << ", e_f " << rec.check.featureError << ", n_- " << rec.check.notCertified
                      << " -> " << (valid ? "ACCEPTED" : "rejected") << " (" << std::fixed << std::setprecision(2)
                      << rec.seconds << " s)";
                    log(1, s.str());
                }
                if (valid)
                {
                    basis = candidate;
                    coefs.swap(candCoefs);
                    ++res.accepted;
                    ++res.acceptedByStage[rec.stage - 1];
                    accepted = true;
                    grendas = (m_opt.selection == Grendas);
                    eligible = working;
                }
                else
                {
                    ++res.rejected;
                    if (group.size() > 1)
                        working.insert(working.end(), group.begin(), group.end());
                }
                res.records.push_back(rec);
            }
        }
    }
    res.basis = basis;
    res.coefs = coefs;
    res.finalSize = basis.size();
    res.check = check(basis, coefs, true);
    res.seconds = clock.stop();
    res.timing.total = res.seconds;
    {
        std::ostringstream s;
        s << "gsMPBESUnrefinement: " << res.initialSize << " -> " << res.finalSize << " MPBES functions, "
          << res.accepted << " accepted (" << res.acceptedByStage[0] << " LS, " << res.acceptedByStage[1]
          << " LO, " << res.acceptedByStage[2] << " NLO), " << res.rejected << " rejected; final e_g "
          << res.check.globalError << ", e_f " << res.check.featureError << ", "
          << (res.check.certified ? "certified" : "NOT certified") << "; " << std::fixed << std::setprecision(1)
          << res.seconds << " s";
        log(1, s.str());
    }
    return res;
}

} // namespace gismo
