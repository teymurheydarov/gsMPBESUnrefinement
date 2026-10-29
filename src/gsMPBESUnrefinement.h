/** @file gsMPBESUnrefinement.h

    @brief The unrefinement algorithm for planar multipatch
    parameterisations in the MPBES basis (Algorithm 1 of Heydarov,
    Buffa, Jüttler): coarsening of the subdomain hierarchy with
    least-squares fitting, linear and non-linear optimisation, under
    regularity and approximation conditions.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): T. Heydarov
*/

#pragma once

#include <gsCore/gsMultiPatch.h>
#include <gsCore/gsMultiBasis.h>
#include <gsMSplines/gsMappedSpline.h>
#include <gsMPBESUnrefinement/src/gsMPBES.h>
#include <gsMPBESUnrefinement/src/gsMPBESProjection.h>
#include <gsMPBESUnrefinement/src/gsBernsteinCertificate.h>
#include <ostream>

namespace gismo
{

// Forward declaration with the default coefficient type
template<class T = real_t> class gsMPBESUnrefinement;

/**
   @brief Unrefinement of a planar multipatch THB-spline
   parameterisation in the MPBES basis.

   Given a fine parameterisation p0 (THB-spline patches with a
   topology), the algorithm removes cells from the subdomain hierarchy
   level by level, from the finest to the coarsest, and projects p0
   onto every candidate space (least-squares fit, then, if the result
   is irregular, the linear and the non-linear optimisation stage of
   gsMPBESProjection).  A candidate is accepted iff

   - (C1) it is regular: every cell of the hierarchical mesh is
     certified by the Bernstein certificate (gsBernsteinCertificate),
   - (C2) its distance from p0 at the points of a uniform grid on every
     patch is at most epsilonGlobal,
   - (C3) its distance from p0 at equally spaced points on the sides of
     the feature set F is at most epsilonFeature.

   The eligible cells of a level are the cells refined exactly one
   level further (Eq. (4) of the paper).  The level is processed in
   passes over its eligible set, as long as the preceding pass accepted
   a coarsening.  Groups of cells are selected by the geometry-based
   Grendas rule (cells of smallest stretching int_c |J|_F^2 up to a
   share delta of the total) at the start of a pass and after every
   acceptance, and by the lexicographic rule (patch index, then the
   upper right cell first) after a rejection; a rejected group of more
   than one cell returns to the working set.  The fitting subdomain is
   the whole parameterisation (lambda = -1) or the union of the
   supports of the functions meeting the lambda-fold extension of the
   selected cells (Section 6.1), in which case only these functions are
   fitted and the others keep their coefficients.

   The hierarchy of a candidate is produced with
   gsHTensorBasis::unrefineElements, the eligibility test with the
   queries of gsHTree, and the candidate basis is a new gsMPBES.

   \tparam T coefficient type
   \ingroup MSplines
*/
template<class T>
class gsMPBESUnrefinement
{
public:
    enum Selection { Grendas = 0, Lexicographic = 1 };

    /// The sample points of the fit
    enum Samples
    {
        UniformGrid = 0,   ///< the q_p x q_p grid of the paper (and the covering grid of a local fit)
        GaussPoints = 1,   ///< Gauss points of the elements of the hierarchical mesh of the candidate (the default, as the HBJ code)
        GaussLobattoPoints = 2 ///< as GaussPoints with the Gauss-Lobatto rule: the element ends, so the patch sides, are sampled
    };

    struct Options
    {
        T epsilonGlobal  = 1e6;            ///< tolerance of (C2); a huge value leaves it unconstrained
        T epsilonFeature = 1e6;            ///< tolerance of (C3)
        std::vector<patchSide> feature;    ///< the feature set F (patch sides); empty: the whole outer boundary
        bool interfacesInFeature = false;  ///< also measure (C3) on the interfaces between patches
        Selection selection = Grendas;
        T grendasDelta = 0.2;              ///< delta of the Grendas rule
        int lambda = -1;                   ///< locality parameter of the fitting subdomain, -1: global fitting
        int stages = 3;                    ///< 1: least squares only, 2: with the linear, 3: with the non-linear stage
        T sampleFactor = 1;                ///< f in q_p = max(5, ceil(f sqrt(n_p))) points per direction and patch
        Samples samples = GaussPoints;     ///< sample points of the fit; (C2) stays on the uniform grid
        int gaussNodes = 0;                ///< nodes per direction and element for GaussPoints and GaussLobattoPoints, 0: 2 degree + 2 (gsGaussRule(basis, 2, 2) of the HBJ code)
        bool unitWeights = true;           ///< GaussPoints/GaussLobattoPoints: coincident points merged, every point weight 1 (as the HBJ code and the grid); false: the quadrature weights
        index_t featurePoints = 10;        ///< points per side of the feature set
        int certificateDepth = 8;          ///< bisection depth of the Bernstein certificate
        typename gsMPBESProjection<T>::Options projection;
        int verbosity = 1;                 ///< 0 silent, 1 one line per candidate, 2 details
        std::ostream * log = nullptr;      ///< an additional log stream
    };

    /// A cell of level \a level (0-based) of \a patch with lower left
    /// corner (x1, y1) in the index space of that level; eligible cells
    /// are refined to level + 1
    struct Cell
    {
        index_t patch, level, x1, y1;
    };

    /// The acceptance test of a candidate
    struct Check
    {
        T globalError = 0;          ///< e_global
        T featureError = 0;         ///< e_feature
        index_t notCertified = -1;  ///< n_-, -1 if the certificate was not computed
        bool certified = false;
        bool valid = false;         ///< (C1), (C2) and (C3) hold
    };

    /// One candidate
    struct Record
    {
        index_t level = 0;
        index_t pass = 0;           ///< pass over the eligible set of the level (Algorithm 1, line 13)
        std::vector<Cell> group;
        bool viaGrendas = false;
        int stage = 0;              ///< last stage applied: 1 LS, 2 LO, 3 NLO
        bool accepted = false;
        index_t size = 0;           ///< size of the candidate basis
        index_t numFree = 0;        ///< fitted functions
        index_t numSamples = 0;     ///< sample points of the fit
        bool augmented = false;     ///< the grid was rank deficient and the element points were added
        Check check;                ///< of the last stage
        double seconds = 0;         ///< wall-clock time of the candidate
        double secondsBasis = 0;    ///< candidate hierarchy, basis and coefficient transfer
        double secondsAssembly = 0; ///< the fitting system: free set, subdomain, sample points, basis values
        double secondsSolve = 0;    ///< the least-squares solve
        double secondsCheck = 0;    ///< accuracy checks and certificates (of all stages)
        double secondsOptimisation = 0; ///< linear and non-linear stage without their checks
    };

    /// Wall-clock time of the phases, summed over the run
    struct Timing
    {
        double selection = 0;       ///< Grendas groups
        double basis = 0, assembly = 0, solve = 0, check = 0, optimisation = 0;
        double total = 0;
    };

    /// The outcome
    struct Result
    {
        gsMPBES<T>  basis;          ///< the final basis
        gsMatrix<T> coefs;          ///< its control points (size() x 2)
        index_t initialSize = 0, finalSize = 0;
        index_t accepted = 0, rejected = 0;
        index_t acceptedByStage[3] = {0, 0, 0};
        Check   check;              ///< of the final parameterisation
        double  seconds = 0;
        Timing  timing;
        std::vector<Record> records;

        /// The final parameterisation
        gsMappedSpline<2,T> geometry() const { return gsMappedSpline<2,T>(basis, coefs); }
        /// The final parameterisation as THB-spline patches of the coarsened hierarchy
        gsMultiPatch<T> thbPatches() const { return gsMPBESUnrefinement<T>::thbPatches(basis, coefs); }
    };

public:
    /// \param p0 the initial parameterisation (THB-spline patches)
    gsMPBESUnrefinement(const gsMultiPatch<T> & p0, const Options & opt = Options());

    const Options & options() const { return m_opt; }

    /// The feature set as patch sides
    const std::vector<patchSide> & feature() const { return m_feature; }

    /// Runs Algorithm 1
    Result run();

    /// @name Building blocks
    /// @{

    /// The eligible cells of level \a level of all patches, in lexicographic order
    static std::vector<Cell> eligibleCells(const gsMPBES<T> & basis, index_t level);

    /// The hierarchy of \a thb with the cells of \a group removed from the next level
    static gsMultiBasis<T> coarsen(const gsMultiBasis<T> & thb, const std::vector<Cell> & group);

    /// Lexicographic order of the cells: patch, then decreasing v and decreasing u
    static bool lexicographicOrder(const Cell & a, const Cell & b);

    /// The parameter box [u1,u2] x [v1,v2] of a cell
    static std::array<T,4> cellBox(const gsMPBES<T> & basis, const Cell & c);

    /// Control points of a parameterisation of \a from carried to the
    /// basis \a to: a function of \a to whose twin class shares a
    /// component with a function of \a from takes its control point;
    /// the indices of the others (set to zero) go to \a newFunctions
    static gsMatrix<T> transferCoefficients(const gsMPBES<T> & from, const gsMatrix<T> & coefs,
                                            const gsMPBES<T> & to, std::vector<index_t> * newFunctions = nullptr);

    /// A parameterisation in an MPBES basis as THB-spline patches of its hierarchy
    static gsMultiPatch<T> thbPatches(const gsMPBES<T> & basis, const gsMatrix<T> & coefs);

    /// (d+1)^2 interior points of every element of the hierarchical mesh of
    /// every patch that meets \a domain (all elements if \a domain is empty)
    static std::vector< gsMatrix<T> > elementPoints(const gsMPBES<T> & basis,
                                                   const std::vector< std::vector< std::array<T,4> > > & domain);

    /// The acceptance test; the certificate is computed only if
    /// \a certify is set and the accuracy conditions hold
    Check check(const gsMPBES<T> & basis, const gsMatrix<T> & coefs, bool certify = true) const;

    /// The sample points of the fit on \a patch (the grid of (C2))
    const gsMatrix<T> & grid(index_t patch) const { return m_grid[patch]; }

    /// The Grendas group of the working set \a working for the
    /// parameterisation (\a basis, \a coefs)
    std::vector<Cell> grendasSet(const gsMPBES<T> & basis, const gsMatrix<T> & coefs,
                                 const std::vector<Cell> & working) const;
    /// @}

private:
    /// The free functions, the sample points and the fitting domain of a candidate
    void fittingSetup(const gsMPBES<T> & candidate, const std::vector<Cell> & group,
                      const std::vector<index_t> & newFunctions, gsMPBESProjection<T> & proj,
                      index_t & numFree) const;

    /// Gauss (or Gauss-Lobatto) points and weights of the elements of every
    /// patch meeting \a domain (all elements if \a domain is empty)
    void gaussSamples(const gsMPBES<T> & basis, const std::vector< std::vector< std::array<T,4> > > & domain,
                      std::vector< gsMatrix<T> > & uv, std::vector< gsVector<T> > & w) const;

    /// The projection stages on a candidate: true if accepted
    bool project(const gsMPBES<T> & candidate, gsMatrix<T> & coefs, const std::vector<Cell> & group,
                 const std::vector<index_t> & newFunctions, Record & rec) const;

    void log(int level, const std::string & line) const;

private:
    gsMultiPatch<T>         m_p0;
    Options                 m_opt;
    std::vector<patchSide>  m_feature;
    std::vector<index_t>    m_gridSize;       ///< q_p
    std::vector< gsMatrix<T> > m_grid;        ///< per patch, the grid of (C2) and of the fit
    std::vector< gsMatrix<T> > m_gridTarget;  ///< p0 there (n x 2)
    std::vector< gsMatrix<T> > m_featurePts;  ///< per patch, the points of (C3)
    std::vector< gsMatrix<T> > m_featureTarget;
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMPBESUnrefinement.hpp)
#endif
