/** @file gsMPBESUnrefinement_test.cpp

    @brief Tests of the unrefinement algorithm

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): T. Heydarov
*/

#include "gismo_unittest.h"
#include <gsMPBESUnrefinement/src/gsMPBESUnrefinement.h>
#include <gsMPBESUnrefinement/src/gsMPBES.h>
#include <gsMPBESUnrefinement/src/gsMPBESSampler.h>

using namespace gismo;

namespace
{
    // two bilinear patches [0,1]x[0,1] and [1,2]x[0,1] of degree 2, uniformly
    // refined to the given number of levels, as THB-spline patches with topology
    gsMultiPatch<> twoPatches(int levels, bool distort)
    {
        gsMultiPatch<> mp;
        gsKnotVector<> kv(0, 1, 1, 3);   // degree 2, two elements
        gsTensorBSplineBasis<2> tb(kv, kv);
        for (int p = 0; p < 2; ++p)
        {
            gsMatrix<> c(16, 2);
            const real_t g[4] = {0, 1.0 / 3, 2.0 / 3, 1};   // Greville-like, exact bilinear map
            for (int j = 0; j < 4; ++j)
                for (int i = 0; i < 4; ++i)
                    c.row(j * 4 + i) << g[i] + p, g[j];
            gsTHBSplineBasis<2> thb(tb);
            gsTHBSpline<2> patch(thb, c);
            for (int l = 1; l < levels; ++l)   // the whole patch to level l
            {
                const index_t n = 2 << l;
                std::vector<index_t> box = {l, 0, 0, n, n};
                patch.refineElements(box);
            }
            if (distort)
            {
                // the fine control points moved by a smooth map: the geometry is not in
                // the coarser spaces, and the interface stays conforming
                gsMatrix<> & cf = patch.coefs();
                for (index_t i = 0; i != cf.rows(); ++i)
                {
                    const real_t x = cf(i, 0), y = cf(i, 1);
                    cf(i, 0) = x + 0.08 * math::sin(2.5 * y);
                    cf(i, 1) = y + 0.05 * math::cos(3 * x);
                }
            }
            mp.addPatch(patch);
        }
        mp.computeTopology();
        return mp;
    }
}

SUITE(gsMPBESUnrefinement_test)
{
    TEST(eligibleCellsAndCoarsening)
    {
        gsMultiPatch<> mp;
        gsReadFile<>("two_patches_refined.xml", mp);
        gsMPBES<> basis(mp);
        // the level-1 box 4..8 x 0..4 of the left patch: four level-0 cells refined to level 1
        std::vector< gsMPBESUnrefinement<>::Cell > cells = gsMPBESUnrefinement<>::eligibleCells(basis, 0);
        CHECK_EQUAL(4u, cells.size());
        for (const gsMPBESUnrefinement<>::Cell & c : cells)
        {
            CHECK_EQUAL(0, c.patch);
            CHECK_EQUAL(0, c.level);
            CHECK(c.x1 >= 2 && c.x1 <= 3 && c.y1 >= 0 && c.y1 <= 1);
        }
        // lexicographic order: upper right first (decreasing v, then decreasing u)
        CHECK_EQUAL(3, cells[0].x1); CHECK_EQUAL(1, cells[0].y1);
        CHECK_EQUAL(2, cells[1].x1); CHECK_EQUAL(1, cells[1].y1);
        CHECK_EQUAL(3, cells[2].x1); CHECK_EQUAL(0, cells[2].y1);
        CHECK(gsMPBESUnrefinement<>::eligibleCells(basis, 1).empty());
        // removing all four gives the unrefined two-patch basis (66 functions)
        gsMPBES<> coarse(gsMPBESUnrefinement<>::coarsen(basis.thbBases(), cells));
        CHECK_EQUAL(66, coarse.size());
        CHECK_EQUAL(1, coarse.numLevels());
        // removing one cell keeps the other three at level 1
        std::vector< gsMPBESUnrefinement<>::Cell > one(1, cells[0]);
        gsMPBES<> partly(gsMPBESUnrefinement<>::coarsen(basis.thbBases(), one));
        CHECK_EQUAL(3u, gsMPBESUnrefinement<>::eligibleCells(partly, 0).size());
        CHECK(partly.size() < basis.size());
        CHECK(partly.size() > coarse.size());
    }

    TEST(coefficientTransfer)
    {
        gsMultiPatch<> mp;
        gsReadFile<>("two_patches_refined.xml", mp);
        gsMPBES<> basis(mp);
        gsMatrix<> coefs = basis.representation(mp);
        std::vector< gsMPBESUnrefinement<>::Cell > cells = gsMPBESUnrefinement<>::eligibleCells(basis, 0);
        std::vector< gsMPBESUnrefinement<>::Cell > one(1, cells[0]);
        gsMPBES<> cand(gsMPBESUnrefinement<>::coarsen(basis.thbBases(), one));
        std::vector<index_t> fresh;
        gsMatrix<> moved = gsMPBESUnrefinement<>::transferCoefficients(basis, coefs, cand, &fresh);
        CHECK_EQUAL(cand.size(), moved.rows());
        // every level-0 function keeps its control point, new level-0 functions (none here:
        // the level-0 functions were all active) are absent, level-1 functions of the removed
        // cell's area disappear
        index_t carried = 0;
        for (index_t g = 0; g != cand.size(); ++g)
            if (moved.row(g).norm() > 0) ++carried;
        CHECK_EQUAL(cand.size() - static_cast<index_t>(fresh.size()), carried);
        // the carried control points are those of the identical functions
        const gsMatrix<> exact = cand.representation(mp);
        real_t worst = 0;
        for (index_t g = 0; g != cand.size(); ++g)
            if (cand.levelOf(g) == 1)
                worst = math::max(worst, (moved.row(g) - exact.row(g)).norm());
        CHECK(worst < 1e-10);
    }

    TEST(exactGeometryCoarsensCompletely)
    {
        // a bilinear geometry refined to three levels lies in the level-0 space: every
        // candidate is accepted by the least-squares fit and the final basis is level 0
        gsMultiPatch<> mp = twoPatches(3, false);
        gsMPBESUnrefinement<>::Options opt;
        opt.epsilonGlobal = 1e-9;
        opt.epsilonFeature = 1e-9;
        opt.verbosity = 0;
        gsMPBESUnrefinement<> algo(mp, opt);
        CHECK_EQUAL(6u, algo.feature().size());    // the outer boundary: 3 sides of each patch
        gsMPBESUnrefinement<>::Result res = algo.run();
        gsMPBES<> level0(gsMPBESUnrefinement<>::coarsen(res.basis.thbBases(), std::vector<gsMPBESUnrefinement<>::Cell>()));
        CHECK_EQUAL(1, res.basis.numLevels());
        CHECK_EQUAL(28, res.finalSize);            // 2 x 16 - 4 shared
        CHECK_EQUAL(0, res.rejected);
        CHECK_EQUAL(res.accepted, res.acceptedByStage[0]);
        CHECK(res.check.valid);
        CHECK(res.check.globalError < 1e-10);
        CHECK(res.initialSize > res.finalSize);
        // the THB-spline patches reproduce the parameterisation
        gsMultiPatch<> out = res.thbPatches();
        CHECK_EQUAL(2u, out.nPatches());
        gsMatrix<> u = gsMPBESSampler<>::uniformGrid(7, 7);
        for (size_t p = 0; p != 2; ++p)
            CHECK((out.patch(p).eval(u) - mp.patch(p).eval(u)).cwiseAbs().maxCoeff() < 1e-10);
    }

    TEST(gaussPointSamples)
    {
        // the fit at the Gauss points of the candidate's elements, weighted by the
        // quadrature weights: the exact geometry still coarsens completely, with
        // (d+1)^2 points per element and no reserve points
        gsMultiPatch<> mp = twoPatches(3, false);
        gsMPBESUnrefinement<>::Options opt;
        opt.epsilonGlobal = 1e-9;
        opt.epsilonFeature = 1e-9;
        opt.samples = gsMPBESUnrefinement<>::GaussPoints;
        opt.verbosity = 0;
        gsMPBESUnrefinement<>::Result res = gsMPBESUnrefinement<>(mp, opt).run();
        CHECK_EQUAL(28, res.finalSize);
        CHECK_EQUAL(0, res.rejected);
        CHECK(res.check.valid);
        for (const gsMPBESUnrefinement<>::Record & r : res.records)
        {
            CHECK(!r.augmented);
            // 9 points per element of the candidate's mesh on both patches, which have
            // at least their 2 x 2 level-0 cells each
            CHECK(r.numSamples >= 9 * 2 * 4);
            CHECK(r.numSamples % 9 == 0);
        }
        // the distorted geometry with a loose tolerance: a valid coarsening as with the grid
        gsMultiPatch<> md = twoPatches(3, true);
        gsMPBESUnrefinement<>::Options loose = opt;
        loose.epsilonGlobal = loose.epsilonFeature = 0.05;
        gsMPBESUnrefinement<>::Result rd = gsMPBESUnrefinement<>(md, loose).run();
        CHECK(rd.finalSize < rd.initialSize);
        CHECK(rd.check.valid);
    }

    TEST(gaussLobattoSamples)
    {
        // the Gauss-Lobatto rule: as many points per element as the Gauss rule,
        // the exact geometry coarsens completely
        gsMultiPatch<> mp = twoPatches(3, false);
        gsMPBESUnrefinement<>::Options opt;
        opt.epsilonGlobal = 1e-9;
        opt.epsilonFeature = 1e-9;
        opt.samples = gsMPBESUnrefinement<>::GaussLobattoPoints;
        opt.verbosity = 0;
        gsMPBESUnrefinement<>::Result res = gsMPBESUnrefinement<>(mp, opt).run();
        CHECK_EQUAL(28, res.finalSize);
        CHECK_EQUAL(0, res.rejected);
        CHECK(res.check.valid);
        for (const gsMPBESUnrefinement<>::Record & r : res.records)
        {
            CHECK(!r.augmented);
            CHECK(r.numSamples % 9 == 0);
        }
    }

    TEST(toleranceStopsCoarsening)
    {
        // a distorted geometry (not in the coarse spaces): a tight tolerance keeps the
        // fine basis, a loose one coarsens, and a feature tolerance acts on the boundary
        gsMultiPatch<> mp = twoPatches(3, true);
        gsMPBESUnrefinement<>::Options tight;
        tight.epsilonGlobal = 1e-12;
        tight.verbosity = 0;
        gsMPBESUnrefinement<>::Result rt = gsMPBESUnrefinement<>(mp, tight).run();
        CHECK_EQUAL(rt.initialSize, rt.finalSize);
        CHECK_EQUAL(0, rt.accepted);

        gsMPBESUnrefinement<>::Options loose;
        loose.epsilonGlobal = 0.05;
        loose.epsilonFeature = 0.05;
        loose.selection = gsMPBESUnrefinement<>::Lexicographic;
        loose.verbosity = 0;
        gsMPBESUnrefinement<>::Result rl = gsMPBESUnrefinement<>(mp, loose).run();
        CHECK(rl.finalSize < rl.initialSize);
        CHECK(rl.check.valid);
        CHECK(rl.check.globalError <= 0.05);
        CHECK(rl.check.featureError <= 0.05);

        // the same with the Grendas selection and a local fitting subdomain
        gsMPBESUnrefinement<>::Options local = loose;
        local.selection = gsMPBESUnrefinement<>::Grendas;
        local.lambda = 1;
        gsMPBESUnrefinement<>::Result rg = gsMPBESUnrefinement<>(mp, local).run();
        CHECK(rg.finalSize < rg.initialSize);
        CHECK(rg.check.valid);
        bool someLocal = false;
        for (const gsMPBESUnrefinement<>::Record & r : rg.records)
            if (r.numFree < r.size) someLocal = true;
        CHECK(someLocal);
    }

    TEST(sampleCoarsensWithoutRejection)
    {
        // the two-patch sample is a level-0 geometry with a level-1 box: every candidate is
        // exact, including those whose uniform grid is rank deficient (the reserve points
        // of the elements are added)
        gsMultiPatch<> mp;
        gsReadFile<>("two_patches_refined.xml", mp);
        gsMPBESUnrefinement<>::Options opt;
        opt.epsilonGlobal = 1e-9;
        opt.epsilonFeature = 1e-9;
        opt.verbosity = 0;
        gsMPBESUnrefinement<>::Result r = gsMPBESUnrefinement<>(mp, opt).run();
        CHECK_EQUAL(76, r.initialSize);
        CHECK_EQUAL(66, r.finalSize);
        CHECK_EQUAL(0, r.rejected);
        CHECK(r.check.valid);
        bool augmented = false;
        for (const gsMPBESUnrefinement<>::Record & rec : r.records) augmented = augmented || rec.augmented;
        CHECK(augmented);
    }

    TEST(featureSubset)
    {
        gsMultiPatch<> mp = twoPatches(2, true);
        gsMPBESUnrefinement<>::Options opt;
        opt.epsilonFeature = 1e-3;
        opt.feature.push_back(patchSide(0, boundary::south));
        opt.feature.push_back(patchSide(0, boundary::east));   // an interface: skipped
        opt.verbosity = 0;
        gsMPBESUnrefinement<> algo(mp, opt);
        CHECK_EQUAL(1u, algo.feature().size());
        gsMPBESUnrefinement<>::Result r = algo.run();
        CHECK(r.check.featureError <= 1e-3);
        CHECK(r.check.valid);
    }
}
