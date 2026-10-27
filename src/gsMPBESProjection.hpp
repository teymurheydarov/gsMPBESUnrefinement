/** @file gsMPBESProjection.hpp

    @brief Implementation of gsMPBESProjection

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): T. Heydarov
*/

#pragma once

#include <gsMPBESUnrefinement/src/gsMPBESProjection.h>
#include <gsMPBESUnrefinement/src/gsMPBESSampler.h>
#include <gsDomain/gsDomain.h>
#include <gsAssembler/gsGaussRule.h>
#include <gsMatrix/gsSparseSolver.h>
#include <casadi/casadi.hpp>

namespace gismo
{

namespace internal
{

/// An Eigen column-major sparse matrix as a CasADi constant
template<class T>
casadi::DM toCasadi(const gsSparseMatrix<T> & A)
{
    gsSparseMatrix<T> C = A;
    C.makeCompressed();
    std::vector<casadi_int> colind(C.outerIndexPtr(), C.outerIndexPtr() + C.cols() + 1);
    std::vector<casadi_int> row(C.innerIndexPtr(), C.innerIndexPtr() + C.nonZeros());
    std::vector<double> val(C.valuePtr(), C.valuePtr() + C.nonZeros());
    return casadi::DM(casadi::Sparsity(C.rows(), C.cols(), colind, row), val, false);
}

/// A dense matrix as a CasADi constant
template<class T>
casadi::DM toCasadi(const gsMatrix<T> & A)
{
    casadi::DM D = casadi::DM::zeros(A.rows(), A.cols());
    for (index_t j = 0; j != A.cols(); ++j)
        for (index_t i = 0; i != A.rows(); ++i)
            D(i, j) = static_cast<double>(A(i, j));
    return D;
}

/// A CasADi sparse result as an Eigen sparse matrix
template<class T>
gsSparseMatrix<T> fromCasadi(const casadi::DM & D)
{
    std::vector<casadi_int> colind, row;
    D.sparsity().get_ccs(colind, row);
    const std::vector<double> & nz = D.nonzeros();
    std::vector< gsEigen::Triplet<T> > trip;
    trip.reserve(nz.size());
    for (casadi_int j = 0; j != D.size2(); ++j)
        for (casadi_int k = colind[j]; k != colind[j+1]; ++k)
            trip.push_back(gsEigen::Triplet<T>(static_cast<index_t>(row[k]), static_cast<index_t>(j), static_cast<T>(nz[k])));
    gsSparseMatrix<T> S(D.size1(), D.size2());
    S.setFromTriplets(trip.begin(), trip.end());
    return S;
}

/// A dense CasADi result (column vector or matrix) as an Eigen matrix
template<class T>
gsMatrix<T> denseFromCasadi(const casadi::DM & D)
{
    const casadi::DM F = casadi::DM::densify(D);
    const std::vector<double> & v = F.nonzeros();    // column-major
    gsMatrix<T> M(F.size1(), F.size2());
    for (casadi_int j = 0; j != F.size2(); ++j)
        for (casadi_int i = 0; i != F.size1(); ++i)
            M(i, j) = static_cast<T>(v[j * F.size1() + i]);
    return M;
}

/**
   The residual vector r(x) of the objective E + Q as a CasADi function
   of the free coefficients x (nFree x 2): the function returns r and
   its Jacobian with respect to vec(x).  The fitting residuals come
   first, then the integrands of the functionals at the quadrature
   nodes, every entry scaled by the square root of its weight, so that
   E + Q = |r|^2.
*/
template<class T>
casadi::Function residualFunction(const gsSparseMatrix<T> & B, const gsMatrix<T> & P,
                                  const gsMatrix<T> & fixedCoefs, const std::vector<index_t> & freeIdx,
                                  const std::vector< gsSparseMatrix<T> > & Q, const gsVector<T> & qw,
                                  T wU, T wL, T wO, T wS, T wE)
{
    using casadi::MX;
    using casadi::DM;
    using casadi::Slice;
    const index_t N = fixedCoefs.rows();
    const index_t nFree = static_cast<index_t>(freeIdx.size());

    MX x = MX::sym("x", nFree, 2);
    // D = fixed + S x
    gsSparseMatrix<T> S(N, nFree);
    {
        std::vector< gsEigen::Triplet<T> > trip;
        for (index_t k = 0; k != nFree; ++k)
            trip.push_back(gsEigen::Triplet<T>(freeIdx[k], k, T(1)));
        S.setFromTriplets(trip.begin(), trip.end());
    }
    gsMatrix<T> fixed = fixedCoefs;
    for (index_t k = 0; k != nFree; ++k)
        fixed.row(freeIdx[k]).setZero();
    MX D = toCasadi(fixed) + MX::mtimes(toCasadi(S), x);

    std::vector<MX> parts;
    // fitting residuals
    parts.push_back(MX::vec(MX::mtimes(toCasadi(B), D) - toCasadi(P)));

    const bool needQ = (wU > 0 || wL > 0 || wO > 0 || wS > 0 || wE > 0);
    if (needQ)
    {
        std::vector<double> sw(qw.size());
        for (index_t k = 0; k != qw.size(); ++k) sw[k] = std::sqrt(static_cast<double>(qw[k]));
        const DM W = DM::diag(DM(sw));          // sqrt of the quadrature weights
        MX Pu  = MX::mtimes(toCasadi(Q[0]), D);
        MX Pv  = MX::mtimes(toCasadi(Q[1]), D);
        auto scaled = [&](const MX & M, T w) { return MX::mtimes(W, M) * std::sqrt(static_cast<double>(w)); };
        if (wL > 0)   // length: |p_u|^2 + |p_v|^2
        {
            parts.push_back(MX::vec(scaled(Pu, wL)));
            parts.push_back(MX::vec(scaled(Pv, wL)));
        }
        if (wU > 0)   // uniformity: |p_uu|^2 + 2 |p_uv|^2 + |p_vv|^2
        {
            MX Puu = MX::mtimes(toCasadi(Q[2]), D);
            MX Pvv = MX::mtimes(toCasadi(Q[3]), D);
            MX Puv = MX::mtimes(toCasadi(Q[4]), D);
            parts.push_back(MX::vec(scaled(Puu, wU)));
            parts.push_back(MX::vec(scaled(Puv, 2 * wU)));
            parts.push_back(MX::vec(scaled(Pvv, wU)));
        }
        if (wO > 0 || wS > 0 || wE > 0)
        {
            const double eps = 1e-30;
            MX dot = Pu(Slice(), 0) * Pv(Slice(), 0) + Pu(Slice(), 1) * Pv(Slice(), 1);
            MX nu  = Pu(Slice(), 0) * Pu(Slice(), 0) + Pu(Slice(), 1) * Pu(Slice(), 1);
            MX nv  = Pv(Slice(), 0) * Pv(Slice(), 0) + Pv(Slice(), 1) * Pv(Slice(), 1);
            if (wO > 0)   // orthogonality: (p_u . p_v)^2
                parts.push_back(scaled(dot, wO));
            if (wS > 0)   // skewness: (p_u . p_v)^2 / (|p_u|^2 |p_v|^2)
                parts.push_back(scaled(dot / MX::sqrt(nu * nv + eps), wS));
            if (wE > 0)   // eccentricity: (p_u . p_uu / |p_u|^2)^2 + (p_v . p_vv / |p_v|^2)^2
            {
                MX Puu = MX::mtimes(toCasadi(Q[2]), D);
                MX Pvv = MX::mtimes(toCasadi(Q[3]), D);
                MX eu = (Pu(Slice(), 0) * Puu(Slice(), 0) + Pu(Slice(), 1) * Puu(Slice(), 1)) / (nu + eps);
                MX ev = (Pv(Slice(), 0) * Pvv(Slice(), 0) + Pv(Slice(), 1) * Pvv(Slice(), 1)) / (nv + eps);
                parts.push_back(scaled(eu, wE));
                parts.push_back(scaled(ev, wE));
            }
        }
    }
    MX r = MX::vertcat(parts);
    MX J = MX::jacobian(r, x);
    return casadi::Function("residual", {x}, {r, J}, {"x"}, {"r", "J"});
}

} // namespace internal

template<class T>
gsMPBESProjection<T>::gsMPBESProjection(const gsMPBES<T> & basis, const gsMultiPatch<T> & target,
                                        const Options & opt)
: m_basis(basis), m_target(target), m_opt(opt), m_augmented(false), m_wholeDomain(true), m_quadReady(false)
{
    GISMO_ENSURE(static_cast<index_t>(target.nPatches()) == static_cast<index_t>(basis.nPatches()),
                 "gsMPBESProjection: the reference parameterisation has " << target.nPatches()
                 << " patches, the basis " << basis.nPatches() << ".");
    const index_t N = basis.size();
    m_free.resize(N);
    for (index_t g = 0; g != N; ++g) m_free[g] = g;

    // functions with a component on a patch side
    m_boundary.assign(N, false);
    std::map< std::pair<index_t,index_t>, std::set<index_t> > onSide;   // (patch, level) -> indices on a side
    for (index_t g = 0; g != N; ++g)
        for (const typename gsMPBES<T>::Component & c : basis.components(g))
        {
            const std::pair<index_t,index_t> key(c.patch, c.level);
            typename std::map< std::pair<index_t,index_t>, std::set<index_t> >::iterator it = onSide.find(key);
            if (it == onSide.end())
            {
                std::set<index_t> s;
                for (int side = 1; side <= 4; ++side)
                {
                    const gsMatrix<index_t> bnd = basis.tensorLevel(c.patch, c.level).boundary(static_cast<boxSide>(side));
                    for (index_t k = 0; k != bnd.size(); ++k) s.insert(bnd(k, 0));
                }
                it = onSide.insert(std::make_pair(key, s)).first;
            }
            if (it->second.count(c.index)) { m_boundary[g] = true; break; }
        }
}

template<class T>
void gsMPBESProjection<T>::sampleTargets(const std::vector< gsMatrix<T> > & uv,
                                         gsSparseMatrix<T> & B, gsMatrix<T> & P) const
{
    GISMO_ASSERT(uv.size() == m_basis.nPatches(), "gsMPBESProjection: one point set per patch.");
    const index_t N = m_basis.size();
    index_t total = 0;
    for (size_t p = 0; p != uv.size(); ++p) total += uv[p].cols();
    P.resize(total, 2);
    std::vector< gsEigen::Triplet<T> > trip;
    index_t row = 0;
    std::vector< gsSparseMatrix<T> > Bp;
    for (size_t p = 0; p != uv.size(); ++p)
    {
        if (uv[p].cols() == 0) continue;
        gsMPBESSampler<T>::sample(m_basis, static_cast<index_t>(p), uv[p], 0, Bp);
        for (index_t j = 0; j != Bp[0].outerSize(); ++j)
            for (typename gsSparseMatrix<T>::InnerIterator it(Bp[0], j); it; ++it)
                trip.push_back(gsEigen::Triplet<T>(row + it.row(), it.col(), it.value()));
        gsMatrix<T> xy;
        m_target.patch(p).eval_into(uv[p], xy);
        P.middleRows(row, uv[p].cols()) = xy.transpose();
        row += uv[p].cols();
    }
    B.resize(total, N);
    B.setFromTriplets(trip.begin(), trip.end());
}

template<class T>
void gsMPBESProjection<T>::setSamples(const std::vector< gsMatrix<T> > & uv)
{
    sampleTargets(uv, m_B, m_P);
    m_Bplain = m_B;
    m_Pplain = m_P;
    m_augmented = false;
}

template<class T>
void gsMPBESProjection<T>::setSamples(const std::vector< gsMatrix<T> > & uv, const std::vector< gsVector<T> > & weights)
{
    sampleTargets(uv, m_Bplain, m_Pplain);
    GISMO_ASSERT(weights.size() == uv.size(), "gsMPBESProjection: one weight vector per patch.");
    gsVector<T> sw(m_Bplain.rows());
    index_t row = 0;
    for (size_t p = 0; p != uv.size(); ++p)
    {
        GISMO_ASSERT(weights[p].size() == uv[p].cols(), "gsMPBESProjection: one weight per sample point.");
        for (index_t k = 0; k != weights[p].size(); ++k)
            sw[row++] = math::sqrt(weights[p][k]);
    }
    m_B = sw.asDiagonal() * m_Bplain;
    m_P = sw.asDiagonal() * m_Pplain;
    m_augmented = false;
}

template<class T>
void gsMPBESProjection<T>::setAugmentation(const std::vector< gsMatrix<T> > & uv)
{
    sampleTargets(uv, m_B2, m_P2);
    m_reserve = PointSource();
    m_augmented = false;
}

template<class T>
void gsMPBESProjection<T>::setAugmentation(const PointSource & source)
{
    m_reserve = source;
    m_B2.resize(0, 0);
    m_P2.resize(0, 2);
    m_augmented = false;
}

template<class T>
void gsMPBESProjection<T>::setFreeFunctions(const std::vector<bool> & free)
{
    const index_t N = m_basis.size();
    m_free.clear();
    if (free.empty())
    {
        for (index_t g = 0; g != N; ++g) m_free.push_back(g);
        return;
    }
    GISMO_ASSERT(static_cast<index_t>(free.size()) == N, "gsMPBESProjection: one flag per function.");
    for (index_t g = 0; g != N; ++g)
        if (free[g]) m_free.push_back(g);
}

template<class T>
void gsMPBESProjection<T>::setFittingDomain(const std::vector< std::vector< std::array<T,4> > > & boxes)
{
    m_domain = boxes;
    m_wholeDomain = boxes.empty();
    m_quadReady = false;
}

template<class T>
void gsMPBESProjection<T>::prepareQuadrature() const
{
    if (m_quadReady) return;
    const index_t N = m_basis.size();
    gsVector<index_t> numNodes(2);
    numNodes.setConstant(m_opt.quadratureNodes);
    gsGaussRule<T> rule(numNodes);
    std::vector< std::vector< gsEigen::Triplet<T> > > trip(5);
    std::vector<T> weights;
    index_t row = 0;
    gsMatrix<T> nodes;
    gsVector<T> w;
    std::vector< gsSparseMatrix<T> > B;
    for (size_t p = 0; p != m_basis.nPatches(); ++p)
    {
        if (!m_wholeDomain && (p >= m_domain.size() || m_domain[p].empty())) continue;
        // the elements of the hierarchical mesh that meet the fitting domain
        std::vector< gsMatrix<T> > elNodes;
        typename gsDomain<T>::Ptr dom = m_basis.thbBasis(static_cast<index_t>(p)).domain();
        for (auto it = dom->beginAll(); it != dom->endAll(); ++it)
        {
            const gsVector<T> lo = it.lowerCorner(), up = it.upperCorner();
            if (!m_wholeDomain)
            {
                bool meets = false;
                for (const std::array<T,4> & b : m_domain[p])
                    if (math::min(b[1], up(0)) > math::max(b[0], lo(0)) &&
                        math::min(b[3], up(1)) > math::max(b[2], lo(1))) { meets = true; break; }
                if (!meets) continue;
            }
            rule.mapTo(lo, up, nodes, w);
            elNodes.push_back(nodes);
            for (index_t k = 0; k != w.size(); ++k) weights.push_back(w[k]);
        }
        if (elNodes.empty()) continue;
        gsMatrix<T> all(2, 0);
        for (const gsMatrix<T> & e : elNodes)
        {
            const index_t c = all.cols();
            all.conservativeResize(2, c + e.cols());
            all.rightCols(e.cols()) = e;
        }
        gsMPBESSampler<T>::sample(m_basis, static_cast<index_t>(p), all, 2, B);
        // B[1] du, B[2] dv, B[3] duu, B[4] dvv, B[5] duv
        for (int q = 0; q != 5; ++q)
            for (index_t j = 0; j != B[q+1].outerSize(); ++j)
                for (typename gsSparseMatrix<T>::InnerIterator it(B[q+1], j); it; ++it)
                    trip[q].push_back(gsEigen::Triplet<T>(row + it.row(), it.col(), it.value()));
        row += all.cols();
    }
    m_Q.resize(5);
    for (int q = 0; q != 5; ++q)
    {
        m_Q[q].resize(row, N);
        m_Q[q].setFromTriplets(trip[q].begin(), trip[q].end());
    }
    m_qw = gsAsConstVector<T>(weights);
    m_quadReady = true;
}

template<class T>
gsMatrix<T> gsMPBESProjection<T>::rowsOf(const gsMatrix<T> & coefs, const std::vector<index_t> & idx)
{
    gsMatrix<T> x(static_cast<index_t>(idx.size()), coefs.cols());
    for (size_t k = 0; k != idx.size(); ++k)
        x.row(k) = coefs.row(idx[k]);
    return x;
}

template<class T>
void gsMPBESProjection<T>::setRows(const gsMatrix<T> & x, const std::vector<index_t> & idx, gsMatrix<T> & coefs)
{
    for (size_t k = 0; k != idx.size(); ++k)
        coefs.row(idx[k]) = x.row(k);
}

template<class T>
std::vector<index_t> gsMPBESProjection<T>::optimised() const
{
    if (!m_opt.freezeBoundary) return m_free;
    std::vector<index_t> idx;
    for (index_t g : m_free)
        if (!m_boundary[g]) idx.push_back(g);
    return idx;
}

template<class T>
void gsMPBESProjection<T>::leastSquares(gsMatrix<T> & coefs) const
{
    GISMO_ENSURE(coefs.rows() == m_basis.size() && coefs.cols() == 2,
                 "gsMPBESProjection: coefficients of size " << m_basis.size() << " x 2 expected.");
    GISMO_ENSURE(m_B.rows() > 0, "gsMPBESProjection: no sample points (call setSamples).");
    const index_t nFree = static_cast<index_t>(m_free.size());
    m_augmented = false;
    if (nFree == 0) return;
    // columns of the free functions, right-hand side without the fixed part
    gsSparseMatrix<T> S(m_basis.size(), nFree);
    {
        std::vector< gsEigen::Triplet<T> > trip;
        for (index_t k = 0; k != nFree; ++k)
            trip.push_back(gsEigen::Triplet<T>(m_free[k], k, T(1)));
        S.setFromTriplets(trip.begin(), trip.end());
    }
    gsMatrix<T> fixed = coefs;
    for (index_t k = 0; k != nFree; ++k) fixed.row(m_free[k]).setZero();

    // the normal equations; a pivot of the factorisation near zero
    // (relative to the largest) reveals a rank-deficient sample set
    auto solve = [&](const gsSparseMatrix<T> & B, const gsMatrix<T> & P, gsMatrix<T> & x) -> bool
    {
        const gsSparseMatrix<T> BS = B * S;
        const gsMatrix<T> rhs = P - B * fixed;
        gsSparseMatrix<T> AtA = BS.transpose() * BS;
        const gsMatrix<T> Atb = BS.transpose() * rhs;
        typename gsSparseSolver<T>::SimplicialLDLT solver;
        solver.compute(AtA);
        if (solver.info() != gsEigen::Success) return false;
        const gsVector<T> d = solver.vectorD();
        const T dmax = d.cwiseAbs().maxCoeff();
        if (!(dmax > 0) || d.cwiseAbs().minCoeff() <= 1e-12 * dmax) return false;
        x = solver.solve(Atb);
        return solver.info() == gsEigen::Success && x.allFinite();
    };
    gsMatrix<T> x;
    bool ok = solve(m_B, m_P, x);
    if (!ok && m_B2.rows() == 0 && m_reserve)
        sampleTargets(m_reserve(), m_B2, m_P2);
    if (!ok && m_B2.rows() > 0)
    {
        // the reserve points join the samples, for this and the following stages
        m_Ball.resize(m_B.rows() + m_B2.rows(), m_basis.size());
        std::vector< gsEigen::Triplet<T> > trip;
        for (index_t j = 0; j != m_B.outerSize(); ++j)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_B, j); it; ++it)
                trip.push_back(gsEigen::Triplet<T>(it.row(), it.col(), it.value()));
        for (index_t j = 0; j != m_B2.outerSize(); ++j)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_B2, j); it; ++it)
                trip.push_back(gsEigen::Triplet<T>(m_B.rows() + it.row(), it.col(), it.value()));
        m_Ball.setFromTriplets(trip.begin(), trip.end());
        m_Pall.resize(m_P.rows() + m_P2.rows(), 2);
        m_Pall << m_P, m_P2;
        m_Bplainall.resize(m_Bplain.rows() + m_B2.rows(), m_basis.size());
        trip.clear();
        for (index_t j = 0; j != m_Bplain.outerSize(); ++j)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_Bplain, j); it; ++it)
                trip.push_back(gsEigen::Triplet<T>(it.row(), it.col(), it.value()));
        for (index_t j = 0; j != m_B2.outerSize(); ++j)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_B2, j); it; ++it)
                trip.push_back(gsEigen::Triplet<T>(m_Bplain.rows() + it.row(), it.col(), it.value()));
        m_Bplainall.setFromTriplets(trip.begin(), trip.end());
        m_Pplainall.resize(m_Pplain.rows() + m_P2.rows(), 2);
        m_Pplainall << m_Pplain, m_P2;
        m_augmented = true;
        ok = solve(m_Ball, m_Pall, x);
        if (m_opt.verbosity > 1)
            gsInfo << "gsMPBESProjection: the sample grid is rank deficient for the " << nFree
                   << " free functions, " << m_B2.rows() << " reserve points added\n";
    }
    if (!ok)
    {
        // minimal-norm least-squares solution by a sparse QR
        const gsSparseMatrix<T> & B = fitMatrix();
        const gsMatrix<T> & P = fitTargets();
        gsEigen::SparseQR< gsSparseMatrix<T>, gsEigen::COLAMDOrdering<index_t> > qr;
        gsSparseMatrix<T> BSc = B * S;
        BSc.makeCompressed();
        qr.compute(BSc);
        x = qr.solve(P - B * fixed);
        if (m_opt.verbosity > 0)
            gsInfo << "gsMPBESProjection: rank-deficient least-squares fit (rank " << qr.rank()
                   << " of " << nFree << "), minimal-norm solution\n";
    }
    setRows(x, m_free, coefs);
}

template<class T>
T gsMPBESProjection<T>::fittingError(const gsMatrix<T> & coefs) const
{
    if (plainMatrix().rows() == 0) return 0;
    const gsMatrix<T> R = plainMatrix() * coefs - plainTargets();
    return R.rowwise().norm().maxCoeff();
}

template<class T>
T gsMPBESProjection<T>::objective(const gsMatrix<T> & coefs, const Weights & w) const
{
    prepareQuadrature();
    casadi::Function F = internal::residualFunction<T>(fitMatrix(), fitTargets(), coefs, m_free, m_Q, m_qw,
                                                       w.uniformity, w.length, w.orthogonality, w.skewness, w.eccentricity);
    std::vector<casadi::DM> out = F(std::vector<casadi::DM>{internal::toCasadi(rowsOf(coefs, m_free))});
    const std::vector<double> & r = out[0].nonzeros();
    double s = 0;
    for (size_t k = 0; k != r.size(); ++k) s += r[k] * r[k];
    return static_cast<T>(s);
}

template<class T>
gsMatrix<T> gsMPBESProjection<T>::gradient(const gsMatrix<T> & coefs, const Weights & w) const
{
    prepareQuadrature();
    casadi::Function F = internal::residualFunction<T>(fitMatrix(), fitTargets(), coefs, m_free, m_Q, m_qw,
                                                       w.uniformity, w.length, w.orthogonality, w.skewness, w.eccentricity);
    std::vector<casadi::DM> out = F(std::vector<casadi::DM>{internal::toCasadi(rowsOf(coefs, m_free))});
    const gsMatrix<T> r = internal::denseFromCasadi<T>(out[0]);
    const gsSparseMatrix<T> J = internal::fromCasadi<T>(out[1]);
    const gsMatrix<T> g = 2 * (J.transpose() * r);     // 2 nFree x 1, column-major over (x, y)
    const index_t nFree = static_cast<index_t>(m_free.size());
    gsMatrix<T> G(nFree, 2);
    G.col(0) = g.topRows(nFree);
    G.col(1) = g.bottomRows(nFree);
    return G;
}

template<class T>
typename gsMPBESProjection<T>::Stats
gsMPBESProjection<T>::gaussNewton(gsMatrix<T> & coefs, const Weights & w, const std::vector<index_t> & idx,
                                  int maxIter, const Validator & valid) const
{
    Stats st;
    const index_t nFree = static_cast<index_t>(idx.size());
    if (nFree == 0) { st.stop = "no free coefficients"; return st; }
    prepareQuadrature();
    casadi::Function F = internal::residualFunction<T>(fitMatrix(), fitTargets(), coefs, idx, m_Q, m_qw,
                                                       w.uniformity, w.length, w.orthogonality, w.skewness, w.eccentricity);
    gsMatrix<T> x = rowsOf(coefs, idx);
    T previousFit = fittingError(coefs);
    int increases = 0;
    st.stop = "iteration limit";
    for (int iter = 0; iter < maxIter; ++iter)
    {
        std::vector<casadi::DM> out = F(std::vector<casadi::DM>{internal::toCasadi(x)});
        const gsMatrix<T> r = internal::denseFromCasadi<T>(out[0]);
        const gsSparseMatrix<T> J = internal::fromCasadi<T>(out[1]);
        gsSparseMatrix<T> H = J.transpose() * J;
        const gsMatrix<T> g = J.transpose() * r;
        typename gsSparseSolver<T>::SimplicialLDLT solver;
        solver.compute(H);
        gsMatrix<T> step;
        if (solver.info() == gsEigen::Success)
            step = solver.solve(g);
        if (solver.info() != gsEigen::Success || !step.allFinite())
        {
            // regularised system
            gsSparseMatrix<T> I(H.rows(), H.cols());
            I.setIdentity();
            const T mu = 1e-10 * math::max(T(1), H.diagonal().maxCoeff());
            H = H + mu * I;
            solver.compute(H);
            step = solver.solve(g);
            if (solver.info() != gsEigen::Success || !step.allFinite())
            {
                st.stop = "singular normal equations";
                break;
            }
        }
        gsMatrix<T> xNew = x;
        xNew.col(0) -= step.topRows(nFree);
        xNew.col(1) -= step.bottomRows(nFree);
        const T stepSize = step.cwiseAbs().maxCoeff(), scale = math::max(T(1), x.cwiseAbs().maxCoeff());
        x = xNew;
        setRows(x, idx, coefs);
        ++st.iterations;
        if (m_opt.verbosity > 1)
            gsInfo << "  Gauss-Newton " << st.iterations << ": |r|^2 = " << r.squaredNorm()
                   << ", step " << stepSize << ", fitting error " << fittingError(coefs) << "\n";
        if (stepSize <= m_opt.stepTolerance * scale)
        {
            st.stop = "negligible step";
            break;
        }
        if (valid && valid(coefs))
        {
            st.valid = true;
            st.stop = "valid iterate";
            break;
        }
        const T fit = fittingError(coefs);
        increases = (fit > previousFit + 1e-12 * previousFit) ? increases + 1 : 0;
        previousFit = fit;
        if (increases >= 2)
        {
            st.stop = "fitting residual increased twice";
            break;
        }
    }
    st.fittingError = fittingError(coefs);
    return st;
}

template<class T>
typename gsMPBESProjection<T>::Stats
gsMPBESProjection<T>::ipopt(gsMatrix<T> & coefs, const Weights & w, const std::vector<index_t> & idx) const
{
    using casadi::MX;
    using casadi::DM;
    Stats st;
    const index_t nFree = static_cast<index_t>(idx.size());
    if (nFree == 0) { st.stop = "no free coefficients"; return st; }
    prepareQuadrature();
    casadi::Function F = internal::residualFunction<T>(fitMatrix(), fitTargets(), coefs, idx, m_Q, m_qw,
                                                       w.uniformity, w.length, w.orthogonality, w.skewness, w.eccentricity);
    MX x = MX::sym("x", 2 * nFree);                       // nlpsol needs a dense vector variable
    MX r = F(std::vector<MX>{MX::reshape(x, nFree, 2)})[0];
    MX f = MX::dot(r, r);
    casadi::Dict opts;
    opts["ipopt.print_level"] = 0;
    opts["print_time"] = 0;
    opts["ipopt.sb"] = "yes";
    opts["ipopt.max_iter"] = m_opt.ipoptIterations;
    opts["ipopt.hessian_approximation"] = "limited-memory";
    casadi::Function solver = casadi::nlpsol("solver", "ipopt", casadi::MXDict{{"x", x}, {"f", f}}, opts);
    casadi::DMDict arg;
    arg["x0"] = casadi::DM::vec(internal::toCasadi(rowsOf(coefs, idx)));
    casadi::DMDict res = solver(arg);
    gsMatrix<T> xNew = internal::denseFromCasadi<T>(res.at("x"));
    xNew.resize(nFree, 2);
    setRows(xNew, idx, coefs);
    casadi::Dict stats = solver.stats();
    st.iterations = static_cast<int>(stats.at("iter_count"));
    st.stop = "ipopt: " + static_cast<std::string>(stats.at("return_status"));
    st.fittingError = fittingError(coefs);
    return st;
}

template<class T>
typename gsMPBESProjection<T>::Stats
gsMPBESProjection<T>::linearStep(gsMatrix<T> & coefs, const Weights & w, const std::vector<index_t> & idx) const
{
    GISMO_ASSERT(w.orthogonality == 0 && w.skewness == 0 && w.eccentricity == 0,
                 "gsMPBESProjection: the linear step takes the quadratic functionals only.");
    Stats st;
    const index_t nFree = static_cast<index_t>(idx.size());
    if (nFree == 0) { st.stop = "no free coefficients"; return st; }
    const index_t N = coefs.rows();
    // the columns of the free functions
    gsSparseMatrix<T> S(N, nFree);
    {
        std::vector< gsEigen::Triplet<T> > trip;
        for (index_t k = 0; k != nFree; ++k)
            trip.push_back(gsEigen::Triplet<T>(idx[k], k, T(1)));
        S.setFromTriplets(trip.begin(), trip.end());
    }
    // The residual of residualFunction() without the non-linear functionals is
    // r = A (fixed + S x) - c, linear in x: one Gauss-Newton step from x0 solves
    // (A S)^T (A S) step = (A S)^T r(x0). The normal equations are summed over the
    // blocks of A: the fitting matrix and the quadrature matrices, the latter with
    // the rows scaled by sqrt(weight of the functional * quadrature weight).
    gsSparseMatrix<T> H(nFree, nFree);
    gsMatrix<T> g = gsMatrix<T>::Zero(nFree, 2);
    auto add = [&](const gsSparseMatrix<T> & M, const gsVector<T> * rowScale, const gsMatrix<T> * target)
    {
        const gsSparseMatrix<T> MS = M * S;
        gsMatrix<T> r0 = M * coefs;
        if (target) r0 -= *target;
        const gsSparseMatrix<T> A = rowScale ? gsSparseMatrix<T>(rowScale->asDiagonal() * MS) : MS;
        const gsMatrix<T>       r = rowScale ? gsMatrix<T>(rowScale->asDiagonal() * r0) : r0;
        H += gsSparseMatrix<T>(A.transpose() * A);
        g += A.transpose() * r;
    };
    add(fitMatrix(), nullptr, &fitTargets());
    if (w.length > 0 || w.uniformity > 0)
    {
        prepareQuadrature();
        const gsVector<T> sw = m_qw.cwiseSqrt();
        if (w.length > 0)   // length: |p_u|^2 + |p_v|^2
        {
            const gsVector<T> s = sw * math::sqrt(w.length);
            add(m_Q[0], &s, nullptr);
            add(m_Q[1], &s, nullptr);
        }
        if (w.uniformity > 0)   // uniformity: |p_uu|^2 + 2 |p_uv|^2 + |p_vv|^2
        {
            const gsVector<T> s  = sw * math::sqrt(w.uniformity);
            const gsVector<T> s2 = sw * math::sqrt(2 * w.uniformity);
            add(m_Q[2], &s,  nullptr);
            add(m_Q[4], &s2, nullptr);
            add(m_Q[3], &s,  nullptr);
        }
    }
    // the solve of gaussNewton(), for both coordinates at once
    typename gsSparseSolver<T>::SimplicialLDLT solver;
    solver.compute(H);
    gsMatrix<T> step;
    if (solver.info() == gsEigen::Success)
        step = solver.solve(g);
    if (solver.info() != gsEigen::Success || !step.allFinite())
    {
        // regularised system
        gsSparseMatrix<T> I(H.rows(), H.cols());
        I.setIdentity();
        const T mu = 1e-10 * math::max(T(1), H.diagonal().maxCoeff());
        H = H + mu * I;
        solver.compute(H);
        step = solver.solve(g);
        if (solver.info() != gsEigen::Success || !step.allFinite())
        {
            st.stop = "singular normal equations";
            st.fittingError = fittingError(coefs);
            return st;
        }
    }
    gsMatrix<T> x = rowsOf(coefs, idx);
    const T stepSize = step.cwiseAbs().maxCoeff(), scale = math::max(T(1), x.cwiseAbs().maxCoeff());
    x -= step;
    setRows(x, idx, coefs);
    st.iterations = 1;
    st.stop = (stepSize <= m_opt.stepTolerance * scale) ? "negligible step" : "iteration limit";
    st.fittingError = fittingError(coefs);
    if (m_opt.verbosity > 1)
        gsInfo << "  linear step: step " << stepSize << ", fitting error " << st.fittingError << "\n";
    return st;
}

template<class T>
typename gsMPBESProjection<T>::Stats
gsMPBESProjection<T>::linearOptimisation(gsMatrix<T> & coefs, const Weights & wIn) const
{
    Weights w = wIn;
    w.orthogonality = w.skewness = w.eccentricity = 0;   // the linear stage uses the quadratic functionals only
    const std::vector<index_t> idx = optimised();
    return (m_opt.solver == "ipopt") ? ipopt(coefs, w, idx) : linearStep(coefs, w, idx);
}

template<class T>
typename gsMPBESProjection<T>::Stats
gsMPBESProjection<T>::nonlinearOptimisation(gsMatrix<T> & coefs, const Weights & w,
                                            const Validator & valid) const
{
    const std::vector<index_t> idx = optimised();
    if (m_opt.solver == "ipopt")
    {
        Stats st = ipopt(coefs, w, idx);
        if (valid) st.valid = valid(coefs);
        return st;
    }
    return gaussNewton(coefs, w, idx, m_opt.maxIterations, valid);
}

template<class T>
typename gsMPBESProjection<T>::Weights
gsMPBESProjection<T>::adaptiveWeights(T globalError, T featureError, T rho, index_t nMinus,
                                      T b, bool linear)
{
    const T rmax = math::max(globalError, math::max(rho, featureError));
    const T beta = T(1.35) * (T(1) + math::min(T(10), T(nMinus) / T(10)));
    auto weight = [&](T r) { return beta * b * (T(0.25) + T(0.75) * (rmax > 0 ? r / rmax : T(0))); };
    Weights w;
    w.uniformity = weight(rho);
    w.length     = weight(featureError);
    if (!linear)
    {
        w.orthogonality = weight(featureError);
        w.skewness      = weight(rho);
        w.eccentricity  = weight(rho);
    }
    return w;
}

} // namespace gismo
