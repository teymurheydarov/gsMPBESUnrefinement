/** @file gsMPBESProjection.h

    @brief The projection stages of the unrefinement algorithm: local
    least-squares fitting, linear and non-linear optimisation of a
    parameterisation in an MPBES space, with the optimisation
    functionals formulated in CasADi.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): T. Heydarov
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsMultiPatch.h>
#include <gsMPBESUnrefinement/src/gsMPBES.h>
#include <functional>
#include <string>
#include <array>

namespace gismo
{

// Forward declaration with the default coefficient type
template<class T = real_t> class gsMPBESProjection;

/**
   @brief The three projection stages of Section 6 of Heydarov, Buffa,
   Jüttler, "An unrefinement algorithm for planar multipatch B-splines
   with enhanced smoothness": for a candidate MPBES basis, the
   reference parameterisation p0 and a set of free functions (the
   interior set B_int, all others keep their coefficients),

   - leastSquares(): the fit (6) of p0 at the sample points, a linear
     least-squares problem solved with a sparse Cholesky factorisation;
   - linearOptimisation(): the fitting error plus the quadratic length
     and uniformity functionals (7); its minimiser is found by one
     Gauss-Newton step, which is exact for a quadratic objective;
   - nonlinearOptimisation(): the objective (8) with the orthogonality,
     skewness and eccentricity functionals added, minimised by
     Gauss-Newton iterations with the stopping rules of Section 6.4
     (an iterate that satisfies the acceptance conditions, a negligible
     step, or a fitting residual that increased twice in a row).

   In the non-linear stage, the residual vector of the objective is
   assembled as a CasADi matrix
   expression of the free coefficients: the fitting residuals and the
   integrands of the functionals at the quadrature nodes of the fitting
   domain, each scaled by the square root of its weight.  CasADi
   supplies the (sparse) Jacobian of this vector by automatic
   differentiation, so that no derivative of a functional is coded by
   hand; the Gauss-Newton normal equations are solved with Eigen.  With
   the option solver = "ipopt", the same objective is minimised by IPOPT
   through CasADi's nlpsol instead of the Gauss-Newton iteration.

   The weights of the functionals are chosen by adaptiveWeights()
   (Section 6.4.1).

   \tparam T coefficient type
   \ingroup MSplines
*/
template<class T>
class gsMPBESProjection
{
public:
    struct Options
    {
        T   bLinear        = 2e-4;   ///< constant b of the weights of the linear stage
        T   bNonlinear     = 2e-3;   ///< constant b of the weights of the non-linear stage
        int maxIterations  = 10;     ///< Gauss-Newton iterations of the non-linear stage
        T   stepTolerance  = 1e-12;  ///< a step smaller than this (relative to the coefficients) ends the iteration
        int quadratureNodes = 2;     ///< Gauss nodes per direction and element for the functionals
        std::string solver = "gn";   ///< "gn": Gauss-Newton, "ipopt": IPOPT via CasADi's nlpsol
        int ipoptIterations = 100;   ///< iteration limit of IPOPT
        /// The functions with a component on a side of a patch (outer
        /// boundary or interface) keep their coefficients during the
        /// optimisation stages; this is not part of the paper's
        /// formulation, off by default
        bool freezeBoundary = false;
        int verbosity = 0;
    };

    /// Weights of the functionals: omega_u, omega_l, omega_o, omega_s, omega_e
    struct Weights
    {
        T uniformity = 0, length = 0, orthogonality = 0, skewness = 0, eccentricity = 0;
        bool nonlinear() const { return orthogonality > 0 || skewness > 0 || eccentricity > 0; }
    };

    /// Outcome of an optimisation stage (the objective of the result is not
    /// computed here: it costs a second construction of the CasADi function;
    /// call objective(coefs, w) when it is needed)
    struct Stats
    {
        int iterations = 0;
        T   fittingError = 0;     ///< largest Euclidean deviation from p0 at the samples
        bool valid = false;       ///< the validity callback accepted an iterate
        std::string stop;         ///< why the iteration ended
    };

    /// Decides whether an iterate is acceptable (conditions (C1)-(C3))
    typedef std::function<bool(const gsMatrix<T> &)> Validator;

public:
    /**
       \param basis   the candidate MPBES basis
       \param target  the reference parameterisation p0 (THB-spline patches)
    */
    gsMPBESProjection(const gsMPBES<T> & basis, const gsMultiPatch<T> & target,
                      const Options & opt = Options());

    /// The sample points of the fit, per patch a 2 x n_p matrix (n_p =
    /// 0 for a patch the fitting domain does not touch); the values of
    /// p0 at these points are the targets
    void setSamples(const std::vector< gsMatrix<T> > & uv);

    /// The same with a weight per point (per patch a vector of length
    /// n_p): the fitting residual of a point is scaled by the square
    /// root of its weight, so that quadrature points and weights give
    /// the discrete L2 projection; the fitting error stays unweighted
    void setSamples(const std::vector< gsMatrix<T> > & uv, const std::vector< gsVector<T> > & weights);

    /// Supplies point sets on demand, per patch a 2 x m_p matrix
    typedef std::function< std::vector< gsMatrix<T> >() > PointSource;

    /**
       @brief Reserve sample points, per patch a 2 x m_p matrix, that
       the fit adds to the samples of setSamples() when the latter are
       rank deficient for the free functions (a uniform grid can miss
       the functions of a locally refined hierarchy); the points of the
       elements of the hierarchical mesh are a safe choice.
    */
    void setAugmentation(const std::vector< gsMatrix<T> > & uv);

    /// The same, with the points produced only when they are needed
    void setAugmentation(const PointSource & source);

    /// True if the reserve points were added by the last least-squares fit
    bool augmentedSamples() const { return m_augmented; }

    /// The free functions (B_int); every function is free if \a free is empty
    void setFreeFunctions(const std::vector<bool> & free);

    /// The fitting domain Delta(c): per patch the boxes [u1,u2]x[v1,v2]
    /// whose union it is; the functionals are integrated over the
    /// elements of the hierarchical mesh that meet it. With an empty
    /// argument the whole parameterisation is the domain.
    void setFittingDomain(const std::vector< std::vector< std::array<T,4> > > & boxes);

    const Options & options() const { return m_opt; }

    /// Projection[1]: the least-squares fit; the free rows of \a coefs
    /// are replaced, the others are the fixed exterior coefficients
    void leastSquares(gsMatrix<T> & coefs) const;

    /// Projection[2]: minimises E + omega_u Q_u + omega_l Q_l from \a coefs
    Stats linearOptimisation(gsMatrix<T> & coefs, const Weights & w) const;

    /// Projection[3]: minimises (8) from \a coefs; \a valid is called
    /// with every iterate and ends the iteration when it returns true
    Stats nonlinearOptimisation(gsMatrix<T> & coefs, const Weights & w,
                                const Validator & valid = Validator()) const;

    /// Largest Euclidean distance between p0 and the parameterisation
    /// \a coefs at the sample points
    T fittingError(const gsMatrix<T> & coefs) const;

    /// The objective E(d) + Q(d) (with the quadrature of the functionals)
    T objective(const gsMatrix<T> & coefs, const Weights & w) const;

    /// Its gradient with respect to the free coefficients (nFree x 2)
    gsMatrix<T> gradient(const gsMatrix<T> & coefs, const Weights & w) const;

    /// Number of free coefficients
    index_t numFree() const { return static_cast<index_t>(m_free.size()); }

    /// Number of sample points of the fit (the reserve points included if they were added)
    index_t numSamples() const { return m_B.rows() + (m_augmented ? m_B2.rows() : 0); }

    /**
       @brief The weights of Section 6.4.1.

       \param globalError  e_global of the preceding candidate
       \param featureError e_feature of the preceding candidate
       \param rho          100 n_- / n_J
       \param nMinus       n_-, the number of cells not certified
       \param b            the stage constant (bLinear or bNonlinear)
       \param linear       the linear stage uses omega_u and omega_l only
    */
    static Weights adaptiveWeights(T globalError, T featureError, T rho, index_t nMinus,
                                   T b, bool linear);

private:
    /// Builds the quadrature of the functionals (nodes, weights, basis derivatives)
    void prepareQuadrature() const;

    /// The basis at the points \a uv and p0 there
    void sampleTargets(const std::vector< gsMatrix<T> > & uv, gsSparseMatrix<T> & B, gsMatrix<T> & P) const;

    /// The (weighted) sample matrix and targets of the fit, with the reserve points if added
    const gsSparseMatrix<T> & fitMatrix() const { return m_augmented ? m_Ball : m_B; }
    const gsMatrix<T> & fitTargets() const { return m_augmented ? m_Pall : m_P; }
    /// The unweighted ones, for the fitting error
    const gsSparseMatrix<T> & plainMatrix() const { return m_augmented ? m_Bplainall : m_Bplain; }
    const gsMatrix<T> & plainTargets() const { return m_augmented ? m_Pplainall : m_Pplain; }

    /// The rows \a idx of \a coefs as a matrix, and back
    static gsMatrix<T> rowsOf(const gsMatrix<T> & coefs, const std::vector<index_t> & idx);
    static void setRows(const gsMatrix<T> & x, const std::vector<index_t> & idx, gsMatrix<T> & coefs);

    /// The free functions of the optimisation stages: the free set
    /// without the boundary functions if they are frozen
    std::vector<index_t> optimised() const;

    /// The Gauss-Newton iteration on the functions \a idx (the non-linear stage)
    Stats gaussNewton(gsMatrix<T> & coefs, const Weights & w, const std::vector<index_t> & idx,
                      int maxIter, const Validator & valid) const;

    /// The single Gauss-Newton step of the linear stage, without CasADi: its
    /// residual is linear in the coefficients, so the Jacobian is assembled
    /// directly from the fitting and quadrature matrices
    Stats linearStep(gsMatrix<T> & coefs, const Weights & w, const std::vector<index_t> & idx) const;

    /// IPOPT on the objective in the functions \a idx
    Stats ipopt(gsMatrix<T> & coefs, const Weights & w, const std::vector<index_t> & idx) const;

private:
    const gsMPBES<T> &      m_basis;
    const gsMultiPatch<T> & m_target;
    Options                 m_opt;

    gsSparseMatrix<T>       m_B;         ///< values of the basis at the samples (stacked over patches), rows scaled by sqrt(weight)
    gsMatrix<T>             m_P;         ///< p0 at the samples, scaled likewise
    gsSparseMatrix<T>       m_Bplain;    ///< the same without the weights
    gsMatrix<T>             m_Pplain;
    mutable gsSparseMatrix<T> m_Bplainall;
    mutable gsMatrix<T>     m_Pplainall;
    mutable gsSparseMatrix<T> m_B2;      ///< the same at the reserve points
    mutable gsMatrix<T>     m_P2;
    PointSource             m_reserve;   ///< produces the reserve points on demand
    mutable bool            m_augmented; ///< the reserve points were added by the last fit
    mutable gsSparseMatrix<T> m_Ball;    ///< samples and reserve points together
    mutable gsMatrix<T>     m_Pall;
    std::vector<index_t>    m_free;      ///< global indices of the free functions
    std::vector<bool>       m_boundary;  ///< functions with a component on a patch side
    std::vector< std::vector< std::array<T,4> > > m_domain;  ///< the fitting domain (empty: everything)
    bool                    m_wholeDomain;

    // quadrature of the functionals (built on demand)
    mutable bool                           m_quadReady;
    mutable gsVector<T>                    m_qw;          ///< weights of all nodes
    mutable std::vector< gsSparseMatrix<T> > m_Q;         ///< Bu, Bv, Buu, Bvv, Buv at the nodes
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMPBESProjection.hpp)
#endif
