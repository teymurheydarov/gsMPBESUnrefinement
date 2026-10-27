/** @file gsMPBES.h

    @brief Multipatch B-splines with enhanced smoothness (MPBES) on
    THB-spline hierarchies: the C0 multipatch extension of truncated
    hierarchical B-splines with activity and truncation decided on the
    domain manifold.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): T. Heydarov
*/

#pragma once

#include <gsMSplines/gsMappedBasis.h>
#include <gsMSplines/gsMappedSpline.h>
#include <gsHSplines/gsTHBSplineBasis.h>
#include <gsHSplines/gsTHBSpline.h>
#include <gsCore/gsMultiBasis.h>
#include <gsCore/gsMultiPatch.h>

namespace gismo
{

// Forward declaration with the default coefficient type (in G+Smo
// itself this line belongs in gsCore/gsForwardDeclarations.h)
template<class T = real_t> class gsMPBES;

/**
   @brief The MPBES basis of a THB-spline hierarchy on a planar
   multipatch domain.

   The construction follows Section 2.2 of Heydarov, Buffa, Jüttler,
   "An unrefinement algorithm for planar THB-spline parameterizations",
   which specialises Buchegger, Jüttler, Mantzaflaris, "Adaptively
   refined multi-patch B-splines with enhanced smoothness" (2016) to
   the C0 case on THB-spline hierarchies:

   - every patch p carries a THB-spline basis, i.e. a subdomain
     hierarchy \f$\Omega^\ell_p\f$ and the tensor-product B-spline bases
     \f$B^\ell_p\f$ of its levels, all patches having the same degree
     and, on every level, uniform knot vectors with the same number of
     knot spans, so that the knot lines match across the interfaces;
   - on every level the boundary functions of neighbouring patches are
     identified across the interfaces through the topology (twins); the
     transitive closure of this relation groups the functions of a level
     into twin classes, which are the MPBES functions \f$B^\ell_M\f$ of
     that level (interior, interface and vertex functions);
   - the Kraft selection is applied on the manifold: a class of level
     \f$\ell\f$ is active iff on every patch it lives on its support lies
     in \f$\Omega^\ell\f$ and on at least one patch it is not contained
     in \f$\Omega^{\ell+1}\f$; a component that is non-zero on an
     interface side whose neighbour has no level \f$\ell\f$ counts as
     lying outside \f$\Omega^\ell\f$;
   - the truncation is applied on the manifold: the finer-level
     children of an active function are removed iff their whole twin
     class lies in the subdomain of their level.

   The result is a gsMappedBasis whose local bases are the finest-level
   tensor-product bases of the patches and whose mapping matrix holds
   the truncated representation of every MPBES function.  All
   evaluation, derivatives, active functions, Paraview output and the
   assemblers come from gsMappedBasis; a geometry in the MPBES space is
   a gsMappedSpline.

   \tparam T coefficient type

   \ingroup MSplines
*/
template<class T>
class gsMPBES : public gsMappedBasis<2,T>
{
public:
    typedef gsMappedBasis<2,T>            Base;
    typedef memory::shared_ptr<gsMPBES>   Ptr;
    typedef memory::unique_ptr<gsMPBES>   uPtr;

    /// A component of an MPBES function: the tensor-product B-spline
    /// \a index of level \a level on patch \a patch
    struct Component
    {
        index_t patch;
        index_t level;
        index_t index;
    };

public:

    gsMPBES() : Base(), m_numLevels(0) { }

    /// Constructs the MPBES basis of the THB-spline bases \a thb (one
    /// gsTHBSplineBasis<2,T> per patch) with the topology of \a thb
    explicit gsMPBES(const gsMultiBasis<T> & thb);

    /// Constructs the MPBES basis of the THB-spline patches of \a mp;
    /// the topology of \a mp is computed if it is absent
    explicit gsMPBES(const gsMultiPatch<T> & mp);

    GISMO_CLONE_FUNCTION(gsMPBES)

public:

    /// @name The hierarchy
    /// @{

    /// Number of levels of the hierarchy (the maximal number over the patches)
    index_t numLevels() const { return m_numLevels; }

    /// Finest level present on patch \a p
    index_t maxLevel(index_t p) const { return static_cast<index_t>(m_levels[p].size()) - 1; }

    /// The THB-spline bases of the patches, with the topology
    const gsMultiBasis<T> & thbBases() const { return m_thb; }

    /// The THB-spline basis of patch \a p
    const gsTHBSplineBasis<2,T> & thbBasis(index_t p) const;

    /// The tensor-product B-spline basis of level \a l on patch \a p
    const gsTensorBSplineBasis<2,T> & tensorLevel(index_t p, index_t l) const
    { return m_levels[p][l]; }

    /// @}

    /// @name The MPBES functions
    /// @{

    /// Components (twin class) of the global function \a g
    const std::vector<Component> & components(index_t g) const { return m_components[g]; }

    /// Level of the global function \a g
    index_t levelOf(index_t g) const { return m_level[g]; }

    /// True if the representation of \a g was truncated
    bool isTruncated(index_t g) const { return m_truncated[g]; }

    /// Global index of the tensor-product B-spline \a index of level \a l
    /// on patch \a p, or -1 if its twin class is not active
    index_t globalIndex(index_t p, index_t l, index_t index) const;

    /// Number of global functions of level \a l
    index_t numFunctionsOfLevel(index_t l) const;

    /// Number of truncated global functions
    index_t numTruncated() const;

    /// Number of global functions with exactly \a n components:
    /// n = 1 interior, n = 2 interface, n >= 3 vertex functions
    index_t numFunctionsWithComponents(index_t n) const;

    /// @}

    /// @name Geometries in the MPBES space
    /// @{

    /**
       @brief Global coefficients of the multipatch \a mp in the MPBES
       space, in the least-squares sense: every patch of \a mp (a
       THB-spline or a tensor-product B-spline on the finest level of
       the patch) is written in the finest-level tensor basis and the
       mapping matrix is solved for the global coefficients.  The
       relative residual, zero iff \a mp lies in the MPBES space, is
       returned in \a residual if given.
    */
    gsMatrix<T> representation(const gsMultiPatch<T> & mp, T * residual = nullptr) const;

    /// The gsMappedSpline with the coefficients of representation()
    gsMappedSpline<2,T> geometry(const gsMultiPatch<T> & mp, T * residual = nullptr) const
    { return gsMappedSpline<2,T>(*this, representation(mp, residual)); }

    /// @}

    /// @name Checks
    /// @{

    /// Largest deviation of the sum of all functions from 1 on a grid
    /// of \a samples points per direction on every patch
    T partitionOfUnityDefect(index_t samples = 7) const;

    /// Largest jump across the interfaces of a function with random
    /// coefficients (seed \a seed), sampled at \a samples points per
    /// interface; zero up to rounding for a C0 basis
    T interfaceJump(index_t samples = 25, int seed = 0) const;

    /// @}

    std::ostream & print(std::ostream & os) const override;

private:

    void init(const gsMultiBasis<T> & thb);

    /// Tensor-product bases of all levels of every patch
    void buildLevels();

    /// Twin classes of every level from the interfaces of the topology
    void identifyTwins();

    /// Coarsest level of the subdomain hierarchy of patch \a p met by the
    /// support of function \a index of level \a l; -1 if the patch has
    /// no level \a l or the function is non-zero on an interface side
    /// without a twin on the neighbour
    int componentLevel(index_t p, index_t l, index_t index) const;

    /// True if every component of the twin class of (p, l, index) has
    /// componentLevel >= l, i.e. the class lies in \f$\Omega^l\f$
    bool classInside(index_t p, index_t l, index_t index) const;

    /// Kraft selection on the manifold
    void selectFunctions();

    /// Truncation on the manifold and the mapping matrix
    void truncateAndMap();

    index_t flat(index_t p, index_t l, index_t index) const
    { return m_levelOffset[l][p] + index; }

private:
    gsMultiBasis<T>                                     m_thb;          // THB bases + topology
    index_t                                             m_numLevels;
    std::vector< std::vector< gsTensorBSplineBasis<2,T> > > m_levels;   // [patch][level]
    std::vector< std::vector<index_t> >                 m_levelOffset;  // [level][patch] -> offset in the flat index space
    std::vector< std::vector<index_t> >                 m_patchOfFlat;  // [level][flat]  -> patch
    std::vector< std::vector<index_t> >                 m_classId;      // [level][flat]  -> twin class
    std::vector< std::vector< std::vector<index_t> > >  m_classes;      // [level][class] -> flat indices
    std::vector< std::vector<char> >                    m_noTwin;       // [level][flat]
    std::vector< std::vector<Component> >               m_components;   // [global]
    std::vector<index_t>                                m_level;        // [global]
    std::vector<bool>                                   m_truncated;    // [global]
    std::vector< std::vector< std::vector<index_t> > >  m_globalIndex;  // [patch][level][index]
    mutable std::vector< std::vector<signed char> >     m_insideMemo;   // [level][flat]: -1 unknown, 0, 1
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMPBES.hpp)
#endif
