/** @file gsBernsteinCertificate.h

    @brief Certificate of the sign of the Jacobian determinant of a
    planar (multipatch) parameterisation from the Bernstein
    coefficients of det J on the cells of its hierarchical mesh.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): T. Heydarov
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsMultiBasis.h>
#include <gsCore/gsMultiPatch.h>
#include <gsMSplines/gsMappedBasis.h>
#include <functional>
#include <limits>

namespace gismo
{

// Forward declaration with the default coefficient type
template<class T = real_t> class gsBernsteinCertificate;

/**
   @brief Decides the regularity condition (C1) of the unrefinement
   algorithm, det J > 0, by a Bernstein certificate.

   On every cell of a mesh on which the parameterisation of degree
   (d1, d2) is a single polynomial, det J is a polynomial of degree
   (2 d1 - 1, 2 d2 - 1).  Its Bernstein coefficients on the cell are
   obtained by interpolation at the 2 d1 x 2 d2 interior nodes
   ((2i+1)/(4 d1), (2j+1)/(4 d2)) of the cell.  The cell is certified
   if all coefficients are positive (negative on a patch with reversed
   orientation), which guarantees the sign of det J on the closed cell
   (convex hull property).  A cell whose coefficients do not all have
   this sign, while its nodal values do, is bisected in both directions
   and its four parts are tested in the same way, each part evaluated
   afresh, up to the depth maxDepth(); a nodal value of the wrong sign
   fails the cell at once, and a cell that is still undecided at the
   depth limit is not certified.  The orientation of a patch is the
   sign of the integral of det J over it, computed from the Bernstein
   coefficients.

   This is the test of Section 4 of Heydarov, Buffa, Jüttler, "An
   unrefinement algorithm for planar multipatch B-splines with enhanced
   smoothness"; notCertified of the result is the number n_- of cells
   that are not certified, which the adaptive weights of the
   optimisation stages use.

   \tparam T coefficient type
   \ingroup MSplines
*/
template<class T>
class gsBernsteinCertificate
{
public:
    /// Outcome on one patch
    struct PatchResult
    {
        int     orientation = 0;    ///< +1 / -1: sign of the signed area, 0: degenerate
        index_t cells       = 0;    ///< cells tested
        index_t certified   = 0;
        index_t failed      = 0;    ///< a nodal value of the wrong sign (or not finite)
        index_t undecided   = 0;    ///< depth limit reached
        index_t subcells    = 0;    ///< (sub)cells tested in total
        int     maxDepth    = 0;    ///< bisection depth used
        T minCoefficient = std::numeric_limits<T>::max(); ///< over the certified (sub)cells: a lower bound of the oriented det J
        T minSample      = std::numeric_limits<T>::max(); ///< over all nodal values
    };

    /// Outcome of the certificate
    struct Result
    {
        bool    certified    = false;   ///< every cell of every patch certified
        index_t cells        = 0;       ///< cells of all patches
        index_t notCertified = 0;       ///< n_-: failed and undecided cells (all cells of a degenerate patch)
        std::vector<PatchResult> patches;
    };

    /// Evaluates det J of the parameterisation on a patch: (patch,
    /// points 2 x n) -> 1 x n
    typedef std::function<void(index_t, const gsMatrix<T> &, gsMatrix<T> &)> DetEvaluator;

public:
    /// \param maxDepth maximal bisection depth of a cell
    explicit gsBernsteinCertificate(int maxDepth = 8) : m_maxDepth(maxDepth) { }

    int maxDepth() const { return m_maxDepth; }
    void setMaxDepth(int depth) { m_maxDepth = depth; }

    /**
       @brief The certificate of a parameterisation given by a mapped
       basis and its coefficients \a coefs (size() x 2).

       \param cells on every patch, a basis whose elements are the cells:
              the parameterisation must be a single polynomial on each
              (the THB-spline bases of the hierarchy of an MPBES basis)
    */
    Result certify(const gsMappedBasis<2,T> & basis, const gsMatrix<T> & coefs,
                   const gsMultiBasis<T> & cells) const;

    /// The certificate of a multipatch geometry (every patch a
    /// polynomial on the elements of its own basis)
    Result certify(const gsMultiPatch<T> & mp) const;

    /**
       @brief The generic certificate.

       \param cells   the bases whose elements are the cells
       \param degrees per patch, the degrees (d1, d2) of the parameterisation
       \param det     the evaluator of det J
    */
    Result certify(const gsMultiBasis<T> & cells,
                   const std::vector< std::pair<short_t,short_t> > & degrees,
                   const DetEvaluator & det) const;

private:
    /// Bernstein collocation matrix at the interior nodes, and its inverse
    static gsMatrix<T> collocationMatrix(index_t n);

    /// Nodes of the cell [a,b] x [c,d] (2 x n1 n2, u fastest)
    static gsMatrix<T> nodes(T a, T b, T c, T d, index_t n1, index_t n2);

    /// +1 certified, -1 failed, 0 undecided
    int certifyCell(const DetEvaluator & det, index_t patch,
                    T a, T b, T c, T d, index_t n1, index_t n2,
                    const gsMatrix<T> & A1inv, const gsMatrix<T> & A2inv,
                    int orientation, int depth, PatchResult & st,
                    const gsMatrix<T> * samples) const;

private:
    int m_maxDepth;
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsBernsteinCertificate.hpp)
#endif
