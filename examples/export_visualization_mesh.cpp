/** @file export_visualization_mesh.cpp

    @brief The final parameterisation of unrefinement_example (<output>_final.xml,
    THB-spline patches with topology) written as the visualisation files of the
    research code of the single-patch algorithm (generateVisualizationMesh in
    poissonTHB_example.cpp):

      <prefix>.txt                  mesh lines: x y patch direction, blank line between lines
      <prefix>_boxes.txt            patch level u0 v0 u1 v1 ix0 iy0 ix1 iy1
      <prefix>_boxes_phys.txt       patch level x0 y0 x1 y1 x2 y2 x3 y3 (sw se ne nw)
      <prefix>_patch_meta.txt       patch mirrored jac_det center tangent directions
      <prefix>_interface_meta.txt   patchA patchB sideA sideB reversed ... (sides U0/U1/V0/V1)
      <prefix>_mesh_meta.txt        degree

    The line order and the number of points per line are those of the research code
    (the boxes of every patch, u-lines then v-lines, degree samples per finest cell,
    n^2 points on a patch side; then one full line per interface side in the order of
    the interface file), so that figures/viz_main.py draws the files as it draws
    those of the research code. The boxes are the leaves
    of the hierarchy, each with its level and its indices in that level. Level 0 may
    consist of several elements per direction (the research code has one per patch);
    its knots are then assumed uniform in [0,1].

    This file is part of the G+Smo library.
*/

#include <gismo.h>
#include <iomanip>

using namespace gismo;

namespace
{
const char * sideName(const boxSide & s)
{
    switch (s.index())
    {
    case 1: return "U0";  // west,  u = 0
    case 2: return "U1";  // east,  u = 1
    case 3: return "V0";  // south, v = 0
    default: return "V1"; // north, v = 1
    }
}

gsMatrix<> pointOnSide(const boxSide & s, real_t t)
{
    gsMatrix<> uv(2, 1);
    switch (s.index())
    {
    case 1: uv << 0, t; break;
    case 2: uv << 1, t; break;
    case 3: uv << t, 0; break;
    default: uv << t, 1; break;
    }
    return uv;
}

gsVector<> unit(const gsVector<> & v)
{
    const real_t n = v.norm();
    return (n > 1e-12) ? gsVector<>(v / n) : gsVector<>(gsVector<>::Zero(2));
}

struct Box { index_t level, x0, y0, x1, y1; };

// a file that cannot be opened (e.g. a path beyond the 260 characters of Windows) is an error
void openOut(std::ofstream & f, const std::string & name)
{
    f.open(name.c_str());
    GISMO_ENSURE(f.is_open(), "Cannot write " << name);
}
}

int main(int argc, char * argv[])
{
    std::string fn, prefix = "output_mesh_final";
    gsCmdLine cmd("Writes the visualisation files of the research code for a multipatch of THB-spline patches.");
    cmd.addString("f", "file",   "Input multipatch (e.g. <output>_final.xml of unrefinement_example)", fn);
    cmd.addString("o", "output", "Prefix of the output files", prefix);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    gsMultiPatch<> mp;
    gsReadFile<>(fn, mp);
    GISMO_ENSURE(mp.nPatches() > 0, "No multipatch in " << fn);
    if (mp.nInterfaces() == 0 && mp.nPatches() > 1)
        mp.computeTopology();
    const index_t nP = mp.nPatches();

    // the leaf boxes of every patch, in the indices of their level
    std::vector< std::vector<Box> > boxes(nP);
    std::vector<index_t> pmax(nP, 0);
    // elements of level 0 per direction: 1 in the box format of the research code; more for
    // inputs refined from a multi-element level 0 (uniform knots in [0,1] assumed)
    std::vector<index_t> mu(nP, 1), mv(nP, 1);
    for (index_t p = 0; p != nP; ++p)
    {
        const gsTHBSplineBasis<2> * thb = dynamic_cast<const gsTHBSplineBasis<2> *>(&mp.patch(p).basis());
        GISMO_ENSURE(thb, "Patch " << p << " is not a THB-spline patch.");
        mu[p] = thb->tensorLevel(0).knots(0).numElements();
        mv[p] = thb->tensorLevel(0).knots(1).numElements();
        gsMatrix<index_t> b1, b2;
        gsVector<index_t> lvl;
        thb->tree().getBoxesInLevelIndex(b1, b2, lvl);
        for (index_t k = 0; k != lvl.size(); ++k)
        {
            Box b = { lvl[k], b1(k, 0), b1(k, 1), b2(k, 0), b2(k, 1) };
            boxes[p].push_back(b);
            pmax[p] = math::max(pmax[p], lvl[k]);
        }
    }
    const index_t degree = mp.patch(0).basis().degree(0);
    const index_t s = degree;  // samples per finest index cell

    // interfaces: midpoints, tangents, orientation
    struct Iface { index_t a, b; boxSide sa, sb; bool rev; };
    std::vector<Iface> ifaces;
    for (const boundaryInterface & bi : mp.interfaces())
    {
        Iface f = { bi.first().patch, bi.second().patch, bi.first().side(), bi.second().side(), false };
        // reversed if the start of side a meets the end of side b
        const gsMatrix<> a0 = mp.patch(f.a).eval(pointOnSide(f.sa, 0));
        const gsMatrix<> b0 = mp.patch(f.b).eval(pointOnSide(f.sb, 0));
        const gsMatrix<> b1 = mp.patch(f.b).eval(pointOnSide(f.sb, 1));
        f.rev = (a0 - b1).norm() < (a0 - b0).norm();
        ifaces.push_back(f);
    }

    std::ofstream out;
    out << std::setprecision(12);

    // _mesh_meta.txt
    openOut(out, prefix + "_mesh_meta.txt");
    out << "# Mesh metadata for Python visualization\n" << "degree " << degree << "\n";
    out.close();

    // _boxes.txt and _boxes_phys.txt
    std::ofstream phys;
    openOut(phys, prefix + "_boxes_phys.txt");
    phys << std::setprecision(12);
    openOut(out, prefix + "_boxes.txt");
    out << "# Box structure in parameter space [0,1]x[0,1]\n" << "# Format: patch level u0 v0 u1 v1 ix0 iy0 ix1 iy1\n";
    phys << "# Box physical corners\n" << "# Format: patch level x0 y0 x1 y1 x2 y2 x3 y3  (corners: sw se ne nw)\n";
    for (index_t p = 0; p != nP; ++p)
        for (const Box & b : boxes[p])
        {
            const real_t hu = 1.0 / static_cast<real_t>(mu[p] * (index_t(1) << b.level));
            const real_t hv = 1.0 / static_cast<real_t>(mv[p] * (index_t(1) << b.level));
            const real_t u0 = b.x0 * hu, v0 = b.y0 * hv, u1 = b.x1 * hu, v1 = b.y1 * hv;
            out << p << " " << b.level << " " << u0 << " " << v0 << " " << u1 << " " << v1 << " "
                << b.x0 << " " << b.y0 << " " << b.x1 << " " << b.y1 << "\n";
            gsMatrix<> c(2, 4);
            c << u0, u1, u1, u0,
                 v0, v0, v1, v1;
            const gsMatrix<> xy = mp.patch(p).eval(c);
            phys << p << " " << b.level;
            for (index_t k = 0; k != 4; ++k) phys << " " << xy(0, k) << " " << xy(1, k);
            phys << "\n";
        }
    out.close();
    phys.close();

    // _patch_meta.txt: centre, unit Jacobian columns, det J at (0.5, 0.5)
    openOut(out, prefix + "_patch_meta.txt");
    out << "# Patch visualization metadata\n"
        << "# Format: patch mirrored jac_det center_x center_y u_dir_x u_dir_y v_dir_x v_dir_y\n";
    for (index_t p = 0; p != nP; ++p)
    {
        gsMatrix<> uv(2, 1);
        uv << 0.5, 0.5;
        const gsMatrix<> xy = mp.patch(p).eval(uv);
        const gsMatrix<> J = mp.patch(p).jacobian(uv);
        const real_t det = J.determinant();
        const gsVector<> du = unit(J.col(0)), dv = unit(J.col(1));
        out << p << " " << (det < 0 ? 1 : 0) << " " << det << " " << xy(0, 0) << " " << xy(1, 0) << " "
            << du(0) << " " << du(1) << " " << dv(0) << " " << dv(1) << "\n";
    }
    out.close();

    // _interface_meta.txt
    openOut(out, prefix + "_interface_meta.txt");
    out << "# Interface visualization metadata\n"
        << "# Format: patchA patchB sideA sideB reversed mirroredA mirroredB gap xA yA tanAx tanAy xB yB tanBx tanBy\n";
    for (const Iface & f : ifaces)
    {
        const gsMatrix<> ua = pointOnSide(f.sa, 0.5), ub = pointOnSide(f.sb, 0.5);
        const gsMatrix<> xa = mp.patch(f.a).eval(ua), xb = mp.patch(f.b).eval(ub);
        const gsMatrix<> Ja = mp.patch(f.a).jacobian(ua), Jb = mp.patch(f.b).jacobian(ub);
        // the tangent of a u-side is the v-derivative and vice versa
        const gsVector<> ta = unit(Ja.col(f.sa.direction() == 0 ? 1 : 0));
        const gsVector<> tb = unit(Jb.col(f.sb.direction() == 0 ? 1 : 0));
        gsMatrix<> c(2, 1);
        c << 0.5, 0.5;
        const bool ma = mp.patch(f.a).jacobian(c).determinant() < 0;
        const bool mb = mp.patch(f.b).jacobian(c).determinant() < 0;
        out << f.a << " " << f.b << " " << sideName(f.sa) << " " << sideName(f.sb) << " " << (f.rev ? 1 : 0) << " "
            << (ma ? 1 : 0) << " " << (mb ? 1 : 0) << " " << (xa - xb).norm() << " "
            << xa(0, 0) << " " << xa(1, 0) << " " << ta(0) << " " << ta(1) << " "
            << xb(0, 0) << " " << xb(1, 0) << " " << tb(0) << " " << tb(1) << "\n";
    }
    out.close();

    // <prefix>.txt: the lines in the order of generateParametricGrid
    std::vector< gsMatrix<> > uvLines;
    std::vector<index_t> linePatch;
    std::vector<char> lineDir;
    for (index_t p = 0; p != nP; ++p)
    {
        const index_t nfinU = mu[p] * (index_t(1) << pmax[p]), nfinV = mv[p] * (index_t(1) << pmax[p]);
        const real_t hfU = 1.0 / static_cast<real_t>(nfinU), hfV = 1.0 / static_cast<real_t>(nfinV);
        for (const Box & b : boxes[p])
        {
            const index_t sf = index_t(1) << (pmax[p] - b.level);
            const index_t X0 = b.x0 * sf, Y0 = b.y0 * sf, X1 = b.x1 * sf, Y1 = b.y1 * sf;
            for (index_t fy = Y0; fy <= Y1; fy += sf)
            {
                const index_t n0 = (X1 - X0) * s + 1;
                const index_t n = (fy == 0 || fy == nfinV) ? n0 * n0 : n0;
                gsMatrix<> line(2, n);
                for (index_t i = 0; i != n; ++i)
                    line.col(i) << (X0 + (X1 - X0) * static_cast<real_t>(i) / (n - 1)) * hfU, fy * hfV;
                uvLines.push_back(line); linePatch.push_back(p); lineDir.push_back('u');
            }
            for (index_t fx = X0; fx <= X1; fx += sf)
            {
                const index_t n0 = (Y1 - Y0) * s + 1;
                const index_t n = (fx == 0 || fx == nfinU) ? n0 * n0 : n0;
                gsMatrix<> line(2, n);
                for (index_t i = 0; i != n; ++i)
                    line.col(i) << fx * hfU, (Y0 + (Y1 - Y0) * static_cast<real_t>(i) / (n - 1)) * hfV;
                uvLines.push_back(line); linePatch.push_back(p); lineDir.push_back('v');
            }
        }
    }
    // one full line per patch side of an interface (appendAlignedInterfaceLines)
    std::set< std::pair<index_t, int> > emitted;
    for (const Iface & f : ifaces)
    {
        const std::pair<index_t, boxSide> ends[2] = { std::make_pair(f.a, f.sa), std::make_pair(f.b, f.sb) };
        for (const std::pair<index_t, boxSide> & e : ends)
        {
            if (!emitted.insert(std::make_pair(e.first, static_cast<int>(e.second.index()))).second) continue;
            // a u-side runs in v and vice versa
            const index_t n = (e.second.direction() == 0 ? mv[e.first] : mu[e.first]) * (index_t(1) << pmax[e.first]) * s + 1;
            gsMatrix<> line(2, n);
            for (index_t i = 0; i != n; ++i)
                line.col(i) = pointOnSide(e.second, static_cast<real_t>(i) / (n - 1));
            uvLines.push_back(line); linePatch.push_back(e.first);
            lineDir.push_back(e.second.direction() == 0 ? 'v' : 'u');
        }
    }
    openOut(out, prefix + ".txt");
    out << "# Visualization Mesh\n" << "# Format: x y patch_id direction\n"
        << "# Blank line separates different grid lines\n" << "# Total lines: " << uvLines.size() << "\n\n";
    for (size_t k = 0; k != uvLines.size(); ++k)
    {
        const gsMatrix<> xy = mp.patch(linePatch[k]).eval(uvLines[k]);
        for (index_t j = 0; j != xy.cols(); ++j)
            out << xy(0, j) << " " << xy(1, j) << " " << linePatch[k] << " " << lineDir[k] << "\n";
        out << "\n";
    }
    out.close();

    gsInfo << fn << ": " << nP << " patches, " << mp.nInterfaces() << " interfaces, " << uvLines.size()
           << " lines -> " << prefix << ".txt and companions\n";
    return EXIT_SUCCESS;
}
