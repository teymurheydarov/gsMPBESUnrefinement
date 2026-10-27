/** @file gsMPBESSampler.h

    @brief Sparse matrices of the values and derivatives of a mapped
    basis at point sets, and the point sets (grids, side points) used
    by the unrefinement algorithm.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): T. Heydarov
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsBoundary.h>
#include <gsMSplines/gsMappedBasis.h>
#include <gsMSplines/gsWeightMapper.h>

namespace gismo
{

// Forward declaration with the default coefficient type
template<class T = real_t> struct gsMPBESSampler;

/**
   @brief Evaluation of a gsMappedBasis<2,T> at arbitrary point sets
   as sparse matrices, and the point sets of the unrefinement
   algorithm.

   gsMappedBasis evaluates point sets that share the active functions
   of their first point; the matrices built here evaluate every point
   with its own active functions, so that uniform grids, side points
   and quadrature nodes of several elements can be handled in one
   call.  With the coefficient matrix D (size() x 2) of a geometry in
   the basis, out[0] * D are its values at the points, out[1] * D and
   out[2] * D its first derivatives, and so on.

   \tparam T coefficient type
   \ingroup MSplines
*/
template<class T>
struct gsMPBESSampler
{
    /**
       @brief The global functions of \a mb on \a patch at the points
       \a u (2 x n) as sparse matrices of size n x mb.size().

       \param order derivatives up to this order (0, 1 or 2)
       \param out   out[0]: values; out[1], out[2]: derivatives in u
                    and v; out[3], out[4], out[5]: second derivatives
                    in uu, vv and uv (the order of gsBasis::deriv2_into)
    */
    static void sample(const gsMappedBasis<2,T> & mb, index_t patch,
                       const gsMatrix<T> & u, int order,
                       std::vector< gsSparseMatrix<T> > & out)
    {
        GISMO_ASSERT(order >= 0 && order <= 2, "gsMPBESSampler: derivatives up to order 2.");
        const gsBasis<T> & base = mb.getBase(patch);
        index_t offset = 0;
        for (index_t q = 0; q < patch; ++q)
            offset += mb.getBase(q).size();
        const index_t n = u.cols(), N = mb.size();

        gsMatrix<index_t> act;
        base.active_into(u, act);
        std::vector< gsMatrix<T> > ders;
        base.evalAllDers_into(u, order, ders);

        const typename gsWeightMapper<T>::LToGMatrix & M = mb.getMapper().asMatrix();
        static const index_t nComp[3] = {1, 2, 3};
        const size_t nOut = (order == 0 ? 1 : (order == 1 ? 3 : 6));
        std::vector< std::vector< gsEigen::Triplet<T> > > trip(nOut);
        const index_t na = act.rows();
        for (index_t k = 0; k != n; ++k)
            for (index_t i = 0; i != na; ++i)
            {
                const index_t local = offset + act(i,k);
                for (typename gsWeightMapper<T>::LToGMatrix::InnerIterator it(M, local); it; ++it)
                {
                    const index_t g = it.col();
                    const T w = it.value();
                    size_t o = 0;
                    for (int d = 0; d <= order; ++d)
                        for (index_t c = 0; c != nComp[d]; ++c, ++o)
                        {
                            const T v = w * ders[d](i * nComp[d] + c, k);
                            if (v != 0)
                                trip[o].push_back(gsEigen::Triplet<T>(k, g, v));
                        }
                }
            }
        out.resize(nOut);
        for (size_t o = 0; o != nOut; ++o)
        {
            out[o].resize(n, N);
            out[o].setFromTriplets(trip[o].begin(), trip[o].end());
        }
    }

    /// Values of a geometry with coefficients \a coefs (size() x 2) at
    /// the points \a u of \a patch, one column per point
    static gsMatrix<T> evalGeometry(const gsMappedBasis<2,T> & mb, index_t patch,
                                    const gsMatrix<T> & u, const gsMatrix<T> & coefs)
    {
        std::vector< gsSparseMatrix<T> > B;
        sample(mb, patch, u, 0, B);
        return (B[0] * coefs).transpose();
    }

    /// The nu x nv grid of the box [u0,u1] x [v0,v1] with the boundary
    /// included, as a 2 x (nu nv) matrix (u varies fastest)
    static gsMatrix<T> uniformGrid(index_t nu, index_t nv,
                                   T u0 = 0, T u1 = 1, T v0 = 0, T v1 = 1)
    {
        GISMO_ASSERT(nu > 1 && nv > 1, "gsMPBESSampler: a grid needs at least two points per direction.");
        gsMatrix<T> g(2, nu * nv);
        for (index_t j = 0; j != nv; ++j)
            for (index_t i = 0; i != nu; ++i)
            {
                g(0, j * nu + i) = u0 + (u1 - u0) * T(i) / T(nu - 1);
                g(1, j * nu + i) = v0 + (v1 - v0) * T(j) / T(nv - 1);
            }
        return g;
    }

    /// \a n equally spaced points on the side \a side of the unit
    /// square, the end points included
    static gsMatrix<T> sidePoints(const boxSide & side, index_t n)
    {
        GISMO_ASSERT(n > 1, "gsMPBESSampler: at least two points per side.");
        gsMatrix<T> g(2, n);
        const index_t dir = side.direction();
        const T fixed = side.parameter() ? T(1) : T(0);
        for (index_t i = 0; i != n; ++i)
        {
            g(dir, i)     = fixed;
            g(1 - dir, i) = T(i) / T(n - 1);
        }
        return g;
    }
};

} // namespace gismo
