/** @file gsBernsteinCertificate.hpp

    @brief Implementation of gsBernsteinCertificate

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): T. Heydarov
*/

#pragma once

#include <gsMPBESUnrefinement/src/gsBernsteinCertificate.h>
#include <gsMPBESUnrefinement/src/gsMPBESSampler.h>
#include <gsDomain/gsDomain.h>

namespace gismo
{

template<class T>
gsMatrix<T> gsBernsteinCertificate<T>::collocationMatrix(index_t n)
{
    // A(i,k) = B_k^{n-1}(tau_i) at the interior nodes tau_i = (2i+1)/(2n)
    gsMatrix<T> A(n, n);
    std::vector<T> binom(n, T(1));
    for (index_t k = 1; k < n; ++k)
        binom[k] = binom[k-1] * T(n - k) / T(k);
    for (index_t i = 0; i < n; ++i)
    {
        const T t = T(2 * i + 1) / T(2 * n);
        for (index_t k = 0; k < n; ++k)
            A(i, k) = binom[k] * math::pow(t, static_cast<int>(k)) * math::pow(T(1) - t, static_cast<int>(n - 1 - k));
    }
    return A;
}

template<class T>
gsMatrix<T> gsBernsteinCertificate<T>::nodes(T a, T b, T c, T d, index_t n1, index_t n2)
{
    gsMatrix<T> u(2, n1 * n2);
    for (index_t j = 0; j < n2; ++j)
        for (index_t i = 0; i < n1; ++i)
        {
            u(0, j * n1 + i) = a + (b - a) * T(2 * i + 1) / T(2 * n1);
            u(1, j * n1 + i) = c + (d - c) * T(2 * j + 1) / T(2 * n2);
        }
    return u;
}

template<class T>
int gsBernsteinCertificate<T>::certifyCell(const DetEvaluator & det, index_t patch,
                                           T a, T b, T c, T d, index_t n1, index_t n2,
                                           const gsMatrix<T> & A1inv, const gsMatrix<T> & A2inv,
                                           int orientation, int depth, PatchResult & st,
                                           const gsMatrix<T> * samples) const
{
    gsMatrix<T> S;
    if (samples)
        S = *samples;
    else
    {
        gsMatrix<T> v;
        det(patch, nodes(a, b, c, d, n1, n2), v);
        S.resize(n1, n2);
        for (index_t j = 0; j < n2; ++j)
            for (index_t i = 0; i < n1; ++i)
                S(i, j) = v(0, j * n1 + i);
    }
    S *= T(orientation);
    const gsMatrix<T> C = A1inv * S * A2inv.transpose();   // oriented Bernstein coefficients

    ++st.subcells;
    st.maxDepth = math::max(st.maxDepth, depth);
    if (!S.allFinite() || !C.allFinite())
        return -1;
    const T minSample = S.minCoeff();
    const T minCoef   = C.minCoeff();
    st.minSample = math::min(st.minSample, minSample);
    if (minCoef > 0)
    {
        st.minCoefficient = math::min(st.minCoefficient, minCoef);
        return +1;
    }
    if (minSample <= 0)
        return -1;
    if (depth >= m_maxDepth)
        return 0;

    // bisection in both directions, every part evaluated afresh
    const T m1 = (a + b) / 2, m2 = (c + d) / 2;
    const T quads[4][4] = { {a, m1, c, m2}, {m1, b, c, m2}, {a, m1, m2, d}, {m1, b, m2, d} };
    int result = +1;
    for (int q = 0; q < 4; ++q)
    {
        const int r = certifyCell(det, patch, quads[q][0], quads[q][1], quads[q][2], quads[q][3],
                                  n1, n2, A1inv, A2inv, orientation, depth + 1, st, nullptr);
        if (r < 0) return -1;
        if (r == 0) result = 0;
    }
    return result;
}

template<class T>
typename gsBernsteinCertificate<T>::Result
gsBernsteinCertificate<T>::certify(const gsMultiBasis<T> & cells,
                                   const std::vector< std::pair<short_t,short_t> > & degrees,
                                   const DetEvaluator & det) const
{
    const index_t nP = static_cast<index_t>(cells.nBases());
    GISMO_ASSERT(static_cast<index_t>(degrees.size()) == nP, "gsBernsteinCertificate: one degree pair per patch.");
    Result res;
    res.patches.resize(nP);
    res.certified = (nP > 0);
    for (index_t p = 0; p < nP; ++p)
    {
        PatchResult & st = res.patches[p];
        const index_t n1 = 2 * degrees[p].first;    // number of Bernstein coefficients of det J in u
        const index_t n2 = 2 * degrees[p].second;
        const gsMatrix<T> A1inv = collocationMatrix(n1).fullPivLu().inverse();
        const gsMatrix<T> A2inv = collocationMatrix(n2).fullPivLu().inverse();

        // the cells of the patch
        std::vector< std::array<T,4> > boxes;
        {
            typename gsDomain<T>::Ptr dom = cells.basis(p).domain();
            for (auto it = dom->beginAll(); it != dom->endAll(); ++it)
            {
                const gsVector<T> lo = it.lowerCorner(), up = it.upperCorner();
                boxes.push_back({lo(0), up(0), lo(1), up(1)});
            }
        }
        st.cells = static_cast<index_t>(boxes.size());

        // nodal values of all cells in one evaluator call
        const index_t nPts = n1 * n2;
        gsMatrix<T> uAll(2, st.cells * nPts), detAll;
        for (index_t ci = 0; ci < st.cells; ++ci)
            uAll.middleCols(ci * nPts, nPts) = nodes(boxes[ci][0], boxes[ci][1], boxes[ci][2], boxes[ci][3], n1, n2);
        det(p, uAll, detAll);
        std::vector< gsMatrix<T> > samples(st.cells);
        for (index_t ci = 0; ci < st.cells; ++ci)
        {
            samples[ci].resize(n1, n2);
            for (index_t j = 0; j < n2; ++j)
                for (index_t i = 0; i < n1; ++i)
                    samples[ci](i, j) = detAll(0, ci * nPts + j * n1 + i);
        }

        // orientation: sign of the integral of det J over the patch
        // (the mean of the Bernstein coefficients is the mean of the polynomial)
        T signedArea = 0;
        for (index_t ci = 0; ci < st.cells; ++ci)
        {
            const gsMatrix<T> C = A1inv * samples[ci] * A2inv.transpose();
            if (C.allFinite())
                signedArea += (boxes[ci][1] - boxes[ci][0]) * (boxes[ci][3] - boxes[ci][2]) * C.sum() / T(nPts);
        }
        st.orientation = (signedArea > 0) ? +1 : ((signedArea < 0) ? -1 : 0);

        if (st.orientation == 0)
            st.failed = st.cells;
        else
            for (index_t ci = 0; ci < st.cells; ++ci)
            {
                const int r = certifyCell(det, p, boxes[ci][0], boxes[ci][1], boxes[ci][2], boxes[ci][3],
                                          n1, n2, A1inv, A2inv, st.orientation, 0, st, &samples[ci]);
                if (r > 0) ++st.certified; else if (r < 0) ++st.failed; else ++st.undecided;
            }

        res.cells += st.cells;
        res.notCertified += st.failed + st.undecided;
        res.certified = res.certified && (st.cells > 0) && (st.certified == st.cells);
    }
    return res;
}

template<class T>
typename gsBernsteinCertificate<T>::Result
gsBernsteinCertificate<T>::certify(const gsMappedBasis<2,T> & basis, const gsMatrix<T> & coefs,
                                   const gsMultiBasis<T> & cells) const
{
    GISMO_ASSERT(coefs.rows() == basis.size() && coefs.cols() == 2,
                 "gsBernsteinCertificate: coefficients of size " << basis.size() << " x 2 expected.");
    std::vector< std::pair<short_t,short_t> > degrees;
    for (size_t p = 0; p != cells.nBases(); ++p)
        degrees.push_back(std::make_pair(cells.basis(p).degree(0), cells.basis(p).degree(1)));
    DetEvaluator det = [&basis, &coefs](index_t patch, const gsMatrix<T> & u, gsMatrix<T> & d)
    {
        std::vector< gsSparseMatrix<T> > B;
        gsMPBESSampler<T>::sample(basis, patch, u, 1, B);
        const gsMatrix<T> Ju = B[1] * coefs, Jv = B[2] * coefs;   // n x 2 each
        d = (Ju.col(0).cwiseProduct(Jv.col(1)) - Ju.col(1).cwiseProduct(Jv.col(0))).transpose();
    };
    return certify(cells, degrees, det);
}

template<class T>
typename gsBernsteinCertificate<T>::Result
gsBernsteinCertificate<T>::certify(const gsMultiPatch<T> & mp) const
{
    gsMultiBasis<T> cells(mp);
    std::vector< std::pair<short_t,short_t> > degrees;
    for (size_t p = 0; p != mp.nPatches(); ++p)
        degrees.push_back(std::make_pair(mp.patch(p).degree(0), mp.patch(p).degree(1)));
    DetEvaluator det = [&mp](index_t patch, const gsMatrix<T> & u, gsMatrix<T> & d)
    {
        gsMatrix<T> der;
        mp.patch(patch).deriv_into(u, der);   // rows: x_u, x_v, y_u, y_v
        d = (der.row(0).cwiseProduct(der.row(3)) - der.row(1).cwiseProduct(der.row(2)));
    };
    return certify(cells, degrees, det);
}

} // namespace gismo
