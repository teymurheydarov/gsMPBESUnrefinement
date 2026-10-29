/** @file unrefinement_example.cpp

    @brief Unrefinement of a planar multipatch THB-spline
    parameterisation in the MPBES basis (Heydarov, Buffa, Jüttler)

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): T. Heydarov
*/

#include <gismo.h>
#include <gsMPBESUnrefinement/src/gsMPBESUnrefinement.h>
#include <fstream>

using namespace gismo;

namespace
{
    // "p:side,p:side,..." with side in {west,east,south,north} or 1..4; an item "p" or
    // "p-q" stands for all boundary sides of the patch p or of the patches p..q
    std::vector<patchSide> parseFeature(const std::string & spec, const gsMultiPatch<> & mp)
    {
        std::vector<patchSide> sides;
        std::stringstream ss(spec);
        std::string item;
        while (std::getline(ss, item, ','))
        {
            if (item.empty()) continue;
            const size_t colon = item.find(':');
            if (colon == std::string::npos)
            {
                const size_t dash = item.find('-');
                const index_t p0 = std::stoi(item.substr(0, dash));
                const index_t p1 = (dash == std::string::npos) ? p0 : std::stoi(item.substr(dash + 1));
                for (const patchSide & bs : mp.boundaries())
                    if (bs.patch >= p0 && bs.patch <= p1) sides.push_back(bs);
                continue;
            }
            const index_t patch = std::stoi(item.substr(0, colon));
            std::string s = item.substr(colon + 1);
            for (char & ch : s) ch = static_cast<char>(std::tolower(ch));
            boxSide side;
            if      (s == "west"  || s == "1" || s == "u0") side = boundary::west;
            else if (s == "east"  || s == "2" || s == "u1") side = boundary::east;
            else if (s == "south" || s == "3" || s == "v0") side = boundary::south;
            else if (s == "north" || s == "4" || s == "v1") side = boundary::north;
            else GISMO_ERROR("--feature: unknown side '" << s << "'");
            sides.push_back(patchSide(patch, side));
        }
        return sides;
    }
}

int main(int argc, char *argv[])
{
    std::string fn = "mask_L3_PatchOnly689_noR.xml";
    real_t epsG = 1e6, epsF = 1e6, delta = 0.2, sampleFactor = 1;
    std::string feature, selection = "grendas", solver = "gn", out = "unrefined", logName, samples = "gauss";
    index_t lambda = -1, stages = 3, verbosity = 1, maxIter = 10, gaussNodes = 0, quadNodes = 2;
    bool interfaces = false, noplot = false, freezeBoundary = false, quadWeights = false;

    gsCmdLine cmd("Unrefinement of a planar multipatch THB-spline parameterisation in the MPBES basis "
                  "(Algorithm 1 of Heydarov, Buffa, Juettler), with the optimisation stages in CasADi.");
    cmd.addString("f", "file",      "Input multipatch of THB-spline patches with topology (XML)", fn);
    cmd.addReal  ("g", "epsilon_g", "Tolerance of the global accuracy condition (C2)", epsG);
    cmd.addReal  ("e", "epsilon_f", "Tolerance of the feature accuracy condition (C3)", epsF);
    cmd.addString("F", "feature",   "Feature set: patch:side,patch:side,... (side: west/east/south/north or 1..4), "
                                    "or p / p-q for all boundary sides of the patches; empty: the whole outer boundary", feature);
    cmd.addSwitch("interfaces",     "Measure (C3) on the interfaces between patches as well", interfaces);
    cmd.addString("s", "selection", "Cell selection: grendas or lex", selection);
    cmd.addReal  ("d", "delta",     "delta of the Grendas selection", delta);
    cmd.addInt   ("l", "lambda",    "Locality parameter of the fitting subdomain (-1: global fitting)", lambda);
    cmd.addReal  ("x", "sample_factor", "f in the q_p = max(5, ceil(f sqrt(n_p))) points per direction of the fit and of (C2)", sampleFactor);
    cmd.addString("p", "samples",   "Sample points of the fit: gauss (Gauss points of the elements of the candidate's hierarchical "
                                    "mesh, as the HBJ code), grid (the uniform grid of the first version of the paper) or lobatto "
                                    "(the Gauss-Lobatto rule, which samples the element ends); (C2) stays on the grid", samples);
    cmd.addInt   ("n", "gauss_nodes", "Gauss or Gauss-Lobatto nodes per direction and element of the fit (0: 2 degree + 2, as the HBJ code)", gaussNodes);
    cmd.addSwitch("quadrature_weights", "gauss/lobatto: weight the points by the quadrature weights instead of weight 1 per point", quadWeights);
    cmd.addInt   ("q", "quad_nodes", "Gauss nodes per direction and element for the functionals of the optimisation stages", quadNodes);
    cmd.addInt   ("t", "stages",    "Projection stages: 1 LS, 2 LS+LO, 3 LS+LO+NLO", stages);
    cmd.addString("S", "solver",    "Optimiser of the stages: gn (Gauss-Newton) or ipopt", solver);
    cmd.addInt   ("i", "iterations","Gauss-Newton iterations of the non-linear stage", maxIter);
    cmd.addSwitch("freeze_boundary", "Keep the boundary functions fixed in the optimisation stages (not part of the paper's formulation)", freezeBoundary);
    cmd.addString("o", "output",    "Prefix of the output files", out);
    cmd.addString("L", "log",       "Log file (default: <output>.log)", logName);
    cmd.addInt   ("v", "verbosity", "0 silent, 1 one line per candidate, 2 the stages", verbosity);
    cmd.addSwitch("noplot",         "No Paraview output", noplot);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    gsMultiPatch<> mp;
    gsReadFile<>(fn, mp);
    GISMO_ENSURE(mp.nPatches() > 0, "No multipatch in " << fn);
    if (mp.nInterfaces() == 0 && mp.nPatches() > 1)
        mp.computeTopology();
    gsInfo << "Read " << gsFileManager::find(fn) << ": " << mp.nPatches() << " patches, "
           << mp.nInterfaces() << " interfaces\n";

    gsMPBESUnrefinement<>::Options opt;
    opt.epsilonGlobal = epsG;
    opt.epsilonFeature = epsF;
    opt.feature = parseFeature(feature, mp);
    opt.interfacesInFeature = interfaces;
    opt.selection = (selection == "lex") ? gsMPBESUnrefinement<>::Lexicographic : gsMPBESUnrefinement<>::Grendas;
    opt.grendasDelta = delta;
    opt.lambda = lambda;
    opt.sampleFactor = sampleFactor;
    opt.samples = (samples == "gauss")   ? gsMPBESUnrefinement<>::GaussPoints
                : (samples == "lobatto") ? gsMPBESUnrefinement<>::GaussLobattoPoints : gsMPBESUnrefinement<>::UniformGrid;
    opt.gaussNodes = gaussNodes;
    opt.unitWeights = !quadWeights;
    opt.projection.quadratureNodes = quadNodes;
    opt.stages = stages;
    opt.projection.solver = solver;
    opt.projection.maxIterations = maxIter;
    opt.projection.freezeBoundary = freezeBoundary;
    opt.projection.verbosity = verbosity;
    opt.verbosity = verbosity;
    std::ofstream logFile((logName.empty() ? out + ".log" : logName).c_str());
    opt.log = &logFile;

    gsMPBESUnrefinement<> algo(mp, opt);
    gsMPBESUnrefinement<>::Result res = algo.run();

    gsInfo << "\nResult: " << res.initialSize << " -> " << res.finalSize << " MPBES functions ("
           << res.accepted << " coarsenings accepted: " << res.acceptedByStage[0] << " by LS, "
           << res.acceptedByStage[1] << " by LO, " << res.acceptedByStage[2] << " by NLO; "
           << res.rejected << " rejected), " << res.seconds << " s\n"
           << "  final errors: global " << res.check.globalError << ", feature " << res.check.featureError
           << "; regularity " << (res.check.certified ? "certified" : "NOT certified") << "\n"
           << "  final hierarchy: " << res.basis.numLevels() << " levels, per level:";
    for (index_t l = 0; l != res.basis.numLevels(); ++l) gsInfo << " " << res.basis.numFunctionsOfLevel(l);
    gsInfo << "\n";
    index_t samplePoints = 0;
    for (const gsMPBESUnrefinement<>::Record & r : res.records) samplePoints += r.numSamples;
    gsInfo << "  sample points of all fits: " << samplePoints << "\n"
           << "  time (s): total " << res.timing.total << ", candidate bases " << res.timing.basis
           << ", assembly of the fitting systems " << res.timing.assembly << ", least-squares solves "
           << res.timing.solve << ", checks and certificates " << res.timing.check << ", optimisation stages "
           << res.timing.optimisation << ", Grendas selection " << res.timing.selection << "\n";

    // the coarsened parameterisation as THB-spline patches, in the format of the input
    gsMultiPatch<> result = res.thbPatches();
    gsFileData<> fd;
    fd << result;
    fd.save(out + "_final.xml");
    gsInfo << "Written " << out << "_final.xml (THB-spline patches of the coarsened hierarchy) and " << out << ".log\n";

    if (!noplot)
    {
        gsWriteParaview(mp, out + "_input", 1000, true);
        gsWriteParaview(result, out + "_final", 1000, true);
        gsInfo << "Paraview: " << out << "_input.pvd, " << out << "_final.pvd\n";
    }
    return 0;
}
