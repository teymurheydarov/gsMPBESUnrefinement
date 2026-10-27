/** @file gsBernsteinCertificate_test.cpp

    @brief Tests of the Bernstein certificate of det J > 0

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): T. Heydarov
*/

#include "gismo_unittest.h"
#include <gsMPBESUnrefinement/src/gsBernsteinCertificate.h>
#include <gsMPBESUnrefinement/src/gsMPBES.h>

using namespace gismo;

namespace
{
    // a bilinear patch with the given corners (tensor order: (0,0), (1,0), (0,1), (1,1))
    gsMultiPatch<> bilinear(real_t x0, real_t y0, real_t x1, real_t y1,
                            real_t x2, real_t y2, real_t x3, real_t y3, int refine = 1)
    {
        gsMatrix<> c(4, 2);
        c << x0, y0, x1, y1, x2, y2, x3, y3;
        gsKnotVector<> kv(0, 1, 0, 2);
        gsTensorBSplineBasis<2> b(kv, kv);
        gsMultiPatch<> mp;
        mp.addPatch(gsTensorBSpline<2>(b, c));
        for (int r = 0; r < refine; ++r)
            mp.uniformRefine();
        return mp;
    }
}

SUITE(gsBernsteinCertificate_test)
{
    TEST(regularSquare)
    {
        gsMultiPatch<> mp = bilinear(0, 0, 1, 0, 0, 1, 1, 1, 2);
        gsBernsteinCertificate<> cert;
        gsBernsteinCertificate<>::Result r = cert.certify(mp);
        CHECK(r.certified);
        CHECK_EQUAL(16, r.cells);
        CHECK_EQUAL(0, r.notCertified);
        CHECK_EQUAL(+1, r.patches[0].orientation);
        CHECK_CLOSE(1.0, r.patches[0].minCoefficient, 1e-12);
        CHECK_EQUAL(0, r.patches[0].maxDepth);
    }

    TEST(mirroredSquare)
    {
        // x -> 1 - x: det J = -1 everywhere, certified with orientation -1
        gsMultiPatch<> mp = bilinear(1, 0, 0, 0, 1, 1, 0, 1, 1);
        gsBernsteinCertificate<>::Result r = gsBernsteinCertificate<>().certify(mp);
        CHECK(r.certified);
        CHECK_EQUAL(-1, r.patches[0].orientation);
        CHECK_EQUAL(4, r.patches[0].certified);
    }

    TEST(foldedQuad)
    {
        // the last two corners swapped: a symmetric bow-tie, signed area zero
        gsMultiPatch<> mp = bilinear(0, 0, 1, 0, 1, 1, 0, 1, 2);
        gsBernsteinCertificate<>::Result r = gsBernsteinCertificate<>().certify(mp);
        CHECK(!r.certified);
        CHECK_EQUAL(0, r.patches[0].orientation);
        CHECK_EQUAL(r.cells, r.patches[0].failed);
        CHECK_EQUAL(r.cells, r.notCertified);
        // an asymmetric fold: positive orientation, cells with negative nodal values
        gsMultiPatch<> mp2 = bilinear(0, 0, 1, 0, 1.5, 1, 0, 1, 2);
        gsBernsteinCertificate<>::Result r2 = gsBernsteinCertificate<>().certify(mp2);
        CHECK(!r2.certified);
        CHECK(r2.patches[0].orientation != 0);
        CHECK(r2.patches[0].failed > 0);
        CHECK_EQUAL(r2.patches[0].failed + r2.patches[0].undecided, r2.notCertified);
        CHECK(r2.patches[0].minSample < 0);
    }

    TEST(nearlyDegenerateSide)
    {
        // p(u,v) = (u^2 - e u, v): det J = 2u - e is negative only for u < e/2.
        // The nodal values of the cell and of its first bisections are positive,
        // the Bernstein coefficient at u = 0 is -e: undecided at a small depth
        // limit, failed once a bisection reaches a node with 2u < e.
        const real_t e = 1e-3;
        gsKnotVector<> kv(0, 1, 0, 3);
        gsTensorBSplineBasis<2> b(kv, kv);
        gsMatrix<> c(9, 2);
        const real_t x[3] = {0, -e / 2, 1 - e}, y[3] = {0, 0.5, 1};
        for (int j = 0; j < 3; ++j)
            for (int i = 0; i < 3; ++i)
                c.row(j * 3 + i) << x[i], y[j];
        gsMultiPatch<> mp;
        mp.addPatch(gsTensorBSpline<2>(b, c));
        gsBernsteinCertificate<>::Result r = gsBernsteinCertificate<>(5).certify(mp);
        CHECK(!r.certified);
        CHECK_EQUAL(0, r.patches[0].failed);
        CHECK_EQUAL(1, r.patches[0].undecided);
        CHECK_EQUAL(5, r.patches[0].maxDepth);
        CHECK(r.patches[0].minSample > 0);
        gsBernsteinCertificate<>::Result r2 = gsBernsteinCertificate<>(12).certify(mp);
        CHECK(!r2.certified);
        CHECK_EQUAL(1, r2.patches[0].failed);
        CHECK_EQUAL(8, r2.patches[0].maxDepth);   // the node 2^-8/8 is the first with 2u < e
        CHECK(r2.patches[0].minSample < 0);
    }

    TEST(certifiedAfterBisection)
    {
        // p(u,v) = (u^3/3 - u^2/2 + (1/4 + delta) u, v): det J = (u - 1/2)^2 + delta > 0,
        // but its degree-5 Bernstein coefficients on [0,1] are not all positive
        const real_t delta = 0.01;
        gsKnotVector<> ku(0, 1, 0, 4), kv(0, 1, 0, 2);
        gsTensorBSplineBasis<2> b(ku, kv);
        gsMatrix<> c(8, 2);
        const real_t x[4] = {0, (0.25 + delta) / 3, 2 * delta / 3, 1.0 / 12 + delta};
        for (int j = 0; j < 2; ++j)
            for (int i = 0; i < 4; ++i)
                c.row(j * 4 + i) << x[i], real_t(j);
        gsMultiPatch<> mp;
        mp.addPatch(gsTensorBSpline<2>(b, c));
        gsBernsteinCertificate<>::Result r = gsBernsteinCertificate<>().certify(mp);
        CHECK(r.certified);
        CHECK(r.patches[0].maxDepth > 0);
        CHECK(r.patches[0].minCoefficient > 0);
        CHECK(r.patches[0].minCoefficient <= delta + 1e-12);
    }

    TEST(mpbesGeometry)
    {
        // the two-patch sample of gsMPBES: the geometry in the MPBES basis is certified,
        // and the certificate agrees with the one of the multipatch itself
        gsMultiPatch<> mp;
        gsReadFile<>("two_patches_refined.xml", mp);
        gsMPBES<> basis(mp);
        gsMatrix<> coefs = basis.representation(mp);
        gsBernsteinCertificate<> cert;
        gsBernsteinCertificate<>::Result r = cert.certify(basis, coefs, basis.thbBases());
        gsBernsteinCertificate<>::Result r0 = cert.certify(mp);
        CHECK(r.certified);
        CHECK(r0.certified);
        CHECK_EQUAL(r0.cells, r.cells);
        CHECK_EQUAL(2u, r.patches.size());
        CHECK_CLOSE(r0.patches[0].minCoefficient, r.patches[0].minCoefficient, 1e-8);
        CHECK_CLOSE(r0.patches[1].minCoefficient, r.patches[1].minCoefficient, 1e-8);
        // a fold: move an interior control point of the first patch far away
        gsMatrix<> bad = coefs;
        index_t interior = -1;
        for (index_t g = 0; g != basis.size() && interior < 0; ++g)
            if (basis.components(g).size() == 1 && basis.components(g)[0].patch == 0 && basis.levelOf(g) == 1)
                interior = g;
        CHECK(interior >= 0);
        bad.row(interior) += gsVector<>::Constant(2, 3.0).transpose();
        gsBernsteinCertificate<>::Result rb = cert.certify(basis, bad, basis.thbBases());
        CHECK(!rb.certified);
        CHECK(rb.notCertified > 0);
    }
}
