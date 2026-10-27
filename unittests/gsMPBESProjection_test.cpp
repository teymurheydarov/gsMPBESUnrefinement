/** @file gsMPBESProjection_test.cpp

    @brief Tests of the projection stages (least squares, linear and
    non-linear optimisation with CasADi)

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): T. Heydarov
*/

#include "gismo_unittest.h"
#include <gsMPBESUnrefinement/src/gsMPBESProjection.h>
#include <gsMPBESUnrefinement/src/gsMPBESSampler.h>
#include <gsMPBESUnrefinement/src/gsMPBES.h>

using namespace gismo;

namespace
{
    // the two-patch sample with its control points moved by a non-affine map,
    // so that the quality functionals and their gradients do not vanish
    gsMultiPatch<> distortedSample()
    {
        gsMultiPatch<> mp;
        gsReadFile<>("two_patches_refined.xml", mp);
        for (size_t p = 0; p != mp.nPatches(); ++p)
        {
            gsMatrix<> & c = mp.patch(p).coefs();
            for (index_t i = 0; i != c.rows(); ++i)
            {
                const real_t x = c(i, 0), y = c(i, 1);
                c(i, 0) = x + 0.1 * math::sin(2 * y);
                c(i, 1) = y + 0.1 * math::cos(3 * x);
            }
        }
        return mp;
    }

    std::vector< gsMatrix<> > grids(size_t nPatches, index_t n)
    {
        std::vector< gsMatrix<> > uv(nPatches);
        for (size_t p = 0; p != nPatches; ++p)
            uv[p] = gsMPBESSampler<>::uniformGrid(n, n);
        return uv;
    }
}

SUITE(gsMPBESProjection_test)
{
    TEST(leastSquaresReproduces)
    {
        gsMultiPatch<> mp;
        gsReadFile<>("two_patches_refined.xml", mp);
        gsMPBES<> basis(mp);
        gsMPBESProjection<> proj(basis, mp);
        proj.setSamples(grids(mp.nPatches(), 12));
        gsMatrix<> coefs;
        coefs.setZero(basis.size(), 2);
        proj.leastSquares(coefs);
        CHECK(proj.fittingError(coefs) < 1e-10);
        const gsMatrix<> exact = basis.representation(mp);
        CHECK((coefs - exact).cwiseAbs().maxCoeff() < 1e-8);

        // half of the functions free, the others fixed at their exact values
        std::vector<bool> free(basis.size(), false);
        for (index_t g = 0; g < basis.size(); g += 2) free[g] = true;
        proj.setFreeFunctions(free);
        gsMatrix<> half = exact;
        for (index_t g = 0; g < basis.size(); g += 2) half.row(g).setZero();
        proj.leastSquares(half);
        CHECK(proj.fittingError(half) < 1e-10);
        CHECK((half - exact).cwiseAbs().maxCoeff() < 1e-8);
    }

    TEST(gradientMatchesFiniteDifferences)
    {
        gsMultiPatch<> mp = distortedSample();
        gsMPBES<> basis(mp);
        gsMPBESProjection<>::Options opt;
        opt.freezeBoundary = false;
        gsMPBESProjection<> proj(basis, mp, opt);
        proj.setSamples(grids(mp.nPatches(), 8));
        gsMatrix<> coefs = basis.representation(mp);
        // a perturbation, so that every functional is away from its minimum
        for (index_t g = 0; g != coefs.rows(); ++g)
            coefs.row(g) += 0.05 * gsVector<>::Random(2).transpose();
        gsMPBESProjection<>::Weights w;
        w.uniformity = 0.3; w.length = 0.2; w.orthogonality = 0.4; w.skewness = 0.5; w.eccentricity = 0.6;
        const gsMatrix<> grad = proj.gradient(coefs, w);
        const real_t h = 1e-6;
        real_t worst = 0;
        for (index_t k = 0; k < basis.size(); k += 7)
            for (index_t d = 0; d < 2; ++d)
            {
                gsMatrix<> cp = coefs, cm = coefs;
                cp(k, d) += h; cm(k, d) -= h;
                const real_t fd = (proj.objective(cp, w) - proj.objective(cm, w)) / (2 * h);
                worst = math::max(worst, math::abs(fd - grad(k, d)) / math::max(real_t(1), math::abs(fd)));
            }
        CHECK(worst < 1e-5);
    }

    TEST(linearStageIsExactMinimiser)
    {
        gsMultiPatch<> mp = distortedSample();
        gsMPBES<> basis(mp);
        gsMPBESProjection<>::Options opt;
        opt.freezeBoundary = false;
        gsMPBESProjection<> proj(basis, mp, opt);
        proj.setSamples(grids(mp.nPatches(), 10));
        gsMatrix<> coefs = basis.representation(mp);
        gsMPBESProjection<>::Weights w;
        w.uniformity = 1e-2; w.length = 1e-2;
        const real_t before = proj.objective(coefs, w);
        gsMPBESProjection<>::Stats st = proj.linearOptimisation(coefs, w);
        CHECK_EQUAL(1, st.iterations);
        CHECK(proj.objective(coefs, w) < before);
        // the objective is quadratic: one Gauss-Newton step reaches the minimum
        const gsMatrix<> grad = proj.gradient(coefs, w);
        CHECK(grad.cwiseAbs().maxCoeff() < 1e-8 * math::max(real_t(1), before));
    }

    TEST(nonlinearStageDecreasesObjective)
    {
        gsMultiPatch<> mp = distortedSample();
        gsMPBES<> basis(mp);
        gsMPBESProjection<>::Options opt;
        opt.freezeBoundary = true;
        gsMPBESProjection<> proj(basis, mp, opt);
        proj.setSamples(grids(mp.nPatches(), 10));
        gsMatrix<> coefs = basis.representation(mp);
        gsMPBESProjection<>::Weights w = gsMPBESProjection<>::adaptiveWeights(0.1, 0.05, 2.0, 3, 2e-3, false);
        CHECK(w.nonlinear());
        const real_t before = proj.objective(coefs, w);
        int calls = 0;
        gsMPBESProjection<>::Stats st = proj.nonlinearOptimisation(coefs, w,
            [&calls](const gsMatrix<>&) { ++calls; return false; });
        CHECK(st.iterations >= 1);
        CHECK_EQUAL(st.iterations, calls);
        CHECK(proj.objective(coefs, w) <= before);
        CHECK(!st.valid);
        // the validator can end the iteration at once
        gsMatrix<> again = basis.representation(mp);
        gsMPBESProjection<>::Stats st2 = proj.nonlinearOptimisation(again, w,
            [](const gsMatrix<>&) { return true; });
        CHECK_EQUAL(1, st2.iterations);
        CHECK(st2.valid);
    }

    TEST(adaptiveWeights)
    {
        // r_max = rho = 2: omega_u = beta b, omega_l = beta b (1/4 + 3/4 * 0.5)
        gsMPBESProjection<>::Weights w = gsMPBESProjection<>::adaptiveWeights(0.5, 1.0, 2.0, 15, 2e-4, true);
        const real_t beta = 1.35 * (1 + 1.5);
        CHECK_CLOSE(beta * 2e-4, w.uniformity, 1e-15);
        CHECK_CLOSE(beta * 2e-4 * (0.25 + 0.75 * 0.5), w.length, 1e-15);
        CHECK_EQUAL(0.0, w.orthogonality);
        CHECK_EQUAL(0.0, w.skewness);
        gsMPBESProjection<>::Weights n = gsMPBESProjection<>::adaptiveWeights(0.5, 1.0, 2.0, 200, 2e-3, false);
        CHECK_CLOSE(1.35 * 11 * 2e-3, n.skewness, 1e-15);
        CHECK_CLOSE(n.length, n.orthogonality, 1e-15);
        CHECK_CLOSE(n.uniformity, n.eccentricity, 1e-15);
    }

    TEST(ipoptSolver)
    {
        gsMultiPatch<> mp = distortedSample();
        gsMPBES<> basis(mp);
        gsMPBESProjection<>::Options opt;
        opt.solver = "ipopt";
        opt.freezeBoundary = true;
        gsMPBESProjection<> proj(basis, mp, opt);
        proj.setSamples(grids(mp.nPatches(), 10));
        gsMatrix<> coefs = basis.representation(mp);
        gsMPBESProjection<>::Weights w;
        w.uniformity = 1e-2; w.length = 1e-2;
        const real_t before = proj.objective(coefs, w);
        proj.linearOptimisation(coefs, w);
        const real_t after = proj.objective(coefs, w);
        CHECK(after < before);
        gsMatrix<> gn = basis.representation(mp);
        gsMPBESProjection<>::Options opt2 = opt;
        opt2.solver = "gn";
        gsMPBESProjection<> proj2(basis, mp, opt2);
        proj2.setSamples(grids(mp.nPatches(), 10));
        proj2.linearOptimisation(gn, w);
        // both reach the minimiser of the same quadratic objective
        CHECK_CLOSE(proj2.objective(gn, w), after, 1e-6 * before);
    }
}
