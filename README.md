# gsMPBESUnrefinement

Implementation, input geometries, run script and logs of the experiments of

> T. Heydarov, A. Buffa, B. Jüttler,
> *An Unrefinement Algorithm for Planar Multipatch B-splines with Enhanced Smoothness*.

The code is a module of [G+Smo](https://github.com/gismo/gismo). Given a planar
multipatch THB-spline parameterization, it removes cells from the subdomain
hierarchy and projects the geometry onto the coarser MPBES space, as long as
the result stays regular (C1) and within the global (C2) and feature (C3)
tolerances. The projection uses up to three stages: least-squares fitting
(LS), linear optimization (LO) and non-linear optimization (NLO); the two
optimization stages are formulated with [CasADi](https://web.casadi.org).

## Contents

| Path | Content |
|---|---|
| `src/` | the module: `gsMPBES` (the MPBES basis), `gsMPBESUnrefinement` (the algorithm), `gsMPBESProjection` (LS, LO, NLO), `gsBernsteinCertificate` (regularity check), `gsMPBESSampler` |
| `examples/unrefinement_example.cpp` | the command-line program used for all experiments |
| `examples/export_visualization_mesh.cpp` | writes the mesh lines and cells of a parameterization as text files, the input of the figure scripts |
| `unittests/` | unit tests of the module (20 tests) |
| `filedata/` | `two_patches_refined.xml`, the input of the unit tests (G+Smo adds this folder to its file search path) |
| `data/geometries/` | the five input parameterizations of the paper and two earlier inputs kept for reference |
| `runs/<run>/` | the logs of the 28 runs of Section 7 and of 6 reference runs on the earlier inputs (see below) |
| `run_all.ps1` | reruns the 28 runs of the paper and compares the results with `runs/` (Windows) |
| `figures/` | the scripts and data that draw the figures of the paper (see below) |

## Requirements

The experiments were computed with the following versions:

| Dependency | Version |
|---|---|
| [G+Smo](https://github.com/gismo/gismo) | commit `b4875c5` (version 26.9.0) |
| [CasADi](https://web.casadi.org) | 3.7.2 (conda-forge) |
| Compiler | MSVC 19.29 (Visual Studio 2019), x64, Release |

## Building

The module is a G+Smo optional module; G+Smo clones it into
`gismo/optional/`. CasADi must provide a `casadi-config.cmake`
(conda-forge, the pip wheel and the release archives do).

```sh
git clone https://github.com/gismo/gismo
git -C gismo checkout b4875c5
cmake -S gismo -B build -DCMAKE_BUILD_TYPE=Release \
  -DGISMO_OPTIONAL="gsMPBESUnrefinement" \
  -DgsMPBESUnrefinement_url=https://github.com/teymurheydarov/gsMPBESUnrefinement \
  -Dcasadi_DIR=<CasADi prefix>/lib/cmake/casadi
cmake --build build --config Release --target unrefinement_example
```

With Visual Studio, add `-G "Visual Studio 16 2019" -A x64` to the first
`cmake` call. The executable is `build/bin/Release/unrefinement_example.exe`;
at run time it needs the CasADi library (`casadi.dll`) on the `PATH`.

Unit tests: add `-DGISMO_BUILD_UNITTESTS=ON` to the first `cmake` call, build
the target `unittests` and run

```sh
build/bin/Release/unittests gsMPBESProjection_test gsMPBESUnrefinement_test gsBernsteinCertificate_test
```

## Reproducing the experiments

On Windows, from the module folder (`gismo/optional/gsMPBESUnrefinement`):

```powershell
powershell -ExecutionPolicy Bypass -File .\run_all.ps1 `
    -Exe <build>\bin\Release\unrefinement_example.exe -CasadiBin <folder with casadi.dll>
```

`-ExecutionPolicy Bypass` applies to this call only; Windows blocks `.ps1`
scripts by default. The script writes every run to `output/<run>/` and prints
the final number of MPBES functions next to the number in `runs/<run>/out.log`.
`-Only run1,run2` selects runs. The times in the logs of the 28 runs add up to
about 90 minutes.

On other systems, and for the 6 reference runs, which the script does not
contain, every run can be started with the command line in the first line of
`runs/<run>/cmd_output.txt` (paths relative to the module folder).

## Reproducing the figures

The figures are drawn from the inputs in `data/geometries/` and the results
`runs/<run>/out_final.xml`; nothing is recomputed. Build the export program and
run the script from the module folder:

```powershell
cmake --build build --config Release --target export_visualization_mesh
py figures\make_figures.py --exe <build>\bin\Release\export_visualization_mesh.exe `
    --path <folder with casadi.dll>
```

The export program, like `unrefinement_example`, needs the CasADi library;
`--path` puts its folder on the `PATH`. Names of figures as arguments select
figures. The script writes `figures/output/<name>.png`, with the file names of
the paper:

| Paper | Figures | Drawn by |
|---|---|---|
| Figure 1 | `mask_bands_original`, `yeti_mp2_thb_approximation_fine_L3_original`, `joystick_approximation_fine_L3_NLO_original`, `rotorRE_original` | `viz_main.py` |
| Figure 2 | `yeti_eg001_ef001_dof239`, `yeti_eg003_ef003_dof157`, `yeti_eg006_ef006_dof122` | `viz_main.py` |
| Figure 3 | `yeti_F_ef010_dof178`, `yeti_F_ef030_dof128`, `yeti_F_ef060_dof117` | `viz_main.py` |
| Figure 4 | `rotor_eg0001_ef0001_dof476`, `rotor_nog_ef0001_dof242` | `viz_main.py` |
| Figure 5 | `mask_bands_ls_dof568`, `mask_bands_lo_dof562`, `mask_bands_nlo_full_dof351` | `viz_main.py` |
| Figure 6 | `tri_wedge_original`, `tri_wedge_ls_dof1994`, `tri_wedge_lo_dof349` | `paper_panel.py` |
| Figure 7 | `joystick_grendas_dof512` | `viz_main.py` |
| Figure 8 | `rotor_subdomains_geometry`, `_lambda0`, `_lambda2`, `_global` | `plot_subdomains.py` |
| Figure 9 | `rotor_lambda0_ef00002_dof356`, `rotor_lambda2_ef00002_dof405`, `rotor_global_ef00002_dof478` | `viz_main.py` |

The number after `_dof` is checked against the last line of `runs/<run>/out.log`.
`export_visualization_mesh` writes the mesh lines and cells of a
parameterization in the text format of the research code of the single-patch
algorithm; `viz_main.py` is the drawing script of that code. The triangle has
several elements on level 0, on which `viz_main.py` does not terminate;
`paper_panel.py` draws its mesh lines directly. Figure 8 uses the files in
`figures/data/`: `first_candidate_<run>.txt`, the cells and the fitting
subdomain of the first candidate, written by the runs `rotor_*_ef00002` with
`--dump_first first_candidate.txt` (the final parameterizations of these runs
are identical to those in `runs/`), and `rotorRE_lines.txt`,
`rotorRE_lines_grid.txt`, patch outlines and grids of `rotorRE.xml` written by
the research code.

With Python 3.9.13, numpy 2.0.2 and matplotlib 3.9.4, without scipy
(`viz_main.py` uses scipy only if it is installed), the script reproduces 18
of the 26 pictures of the paper byte for byte: all except the input pictures
of Yeti, joystick and rotor (Figure 1) and the Yeti results of Figures 2 and 3
except `yeti_F_ef010_dof178`. The published versions of these eight were drawn
from the export and the runs of the research code; the script draws them from
the inputs and runs of this repository, with the same final numbers of MPBES
functions and the same cells. Two differences are visible: the patch numbers
are placed slightly differently, and the published rotor input shows patch 22
refined uniformly, while the redrawn one shows its hierarchy in `rotorRE.xml`
(one level above the other patches, two in the box (0.25, 0.75)^2), as stated
in the caption of Figure 8. Other versions of the Python packages may change
the pixels.

## The runs

Every run folder contains `cmd_output.txt` (the command line and the program
output, including the timing of the phases), `out.log` (one line per
candidate coarsening: the stages invoked, the errors, the decision) and
`out_final.xml` (the coarsened parameterization as THB-spline patches).
Absolute paths in `cmd_output.txt` were replaced by paths relative to the
module folder; otherwise the files are as written by the program.

The logs were produced with the code of this repository. It contains two
implementation changes made after the runs behind the numbers of the paper:
the linear stage assembles its Jacobian, which is constant, directly from the
fitting and quadrature matrices instead of differentiating it with CasADi, and
the stages no longer evaluate the objective a second time for a value the
algorithm does not use. Neither changes a result: rerunning the earlier and
the present code gives identical `out.log` files (every candidate, every
decision, every printed error) and final control points that agree to
1.5e-11; only the times differ.

Commit `9bd6a8a` changed the quantity rho of the adaptive weights of the
optimization stages from 100 n_- / n_J (n_J the points of a uniform
regularity grid) to 100 n_- / n_c (n_c the cells of the hierarchical meshes
of all patches), as in the paper. The runs that reach the linear or
non-linear stage (`mask_*`, `triangle_*`) were computed with the new rho;
the logs of the earlier rho are in the history. The other runs are decided
by the least-squares stage alone, where rho is not computed; rerunning them
with the present code gives the same candidates, decisions and errors.

The fits sample the Gauss points of the elements of the candidate's
hierarchical mesh (2 p + 2 per direction and element, every point with
weight 1), as the code of the single-patch algorithm (Heydarov, Buffa,
Jüttler, *An Unrefinement Algorithm for Planar THB-spline Parameterizations*)
does; conditions
(C2) and (C3) are checked on a uniform grid. The first version of this
repository (commit `78e2643`) fitted on the uniform grid as well; that
sampling remains available with `-p grid`, and the logs of that version are
in the history.

| Run | Input | Paper | MPBES functions |
|---|---|---|---|
| `yeti_eg001_ef001`, `yeti_eg003_ef003`, `yeti_eg006_ef006` | Yeti | Table 2, row 1; Figure 2 | 239, 157, 122 |
| `yeti_nof_ef010`, `yeti_nof_ef030`, `yeti_nof_ef060` | Yeti | Table 2, row 2 | 239, 157, 122 |
| `yeti_F_ef010`, `yeti_F_ef030`, `yeti_F_ef060` | Yeti | Table 2, row 3; Figure 3 | 178, 128, 117 |
| `rotor_eg0001_ef0001`, `rotor_nog_ef0001` | rotor | Section 7.1; Figure 4 | 476, 242 |
| `mask_bands_ls_only`, `mask_bands_lo_only`, `mask_bands_nlo_full` | mask | Table 3; Figure 5 | 568, 562, 351 |
| `triangle_wedge_ls_only`, `triangle_wedge_lo_only`, `triangle_wedge_nlo_full` | triangle | Table 3; Figure 6 | 1994, 349, 349 |
| `joystick_grendas`, `joystick_lex` | joystick | Table 4; Figure 7 | 512, 512 |
| `joystick_global`, `joystick_lambda0`, `joystick_lambda1` | joystick | Section 7.2.3 | 512, 503, 503 |
| `rotor_lambda0_ef00002`, `rotor_lambda2_ef00002`, `rotor_global_ef00002` | rotor | Table 5; Figures 8, 9 | 356, 405, 478 |
| `rotor_global`, `rotor_lambda0`, `rotor_lambda2` | rotor | Section 7.2.3 (text) | 168, 168, 168 |

Notes on the settings:

- An unconstrained tolerance is the program's default `1e6`. The runs
  `yeti_nof_*` and `yeti_F_*` use `-g 1.0`, which never binds on the Yeti
  (the largest global error of a regular candidate in these runs is 0.079).
- The mask runs use `-g 1.5 -e 0.65` (0.65 is less than 1 % of the width of the domain) and the
  feature set `-F` of the 18 true boundary sides (outer boundary and the two holes); the file's
  boundary list contains the interface sides as well.
- The mask input `mask_bands_L3.xml` is regular (every cell certified): the outer boundary and the
  interfaces are those of `mask_L3_PatchOnly689_noR.xml` (the earlier input, whose six folded cells
  at the left hole were repaired by moving interior control points), and a compressed band was added
  inside every patch (det J down to 3 % of the patch median). The earlier input and its runs
  `mask_folded_*` (1089, 1089, 225 at `-g 1.2 -e 0.65`, p0 not certified) are kept for reference.
- The triangle (`triangle_wedge_L5.xml`) is a generated input: three bilinear
  patches forming a thin triangle with two wavy sides, 3169 MPBES functions over
  six levels, every cell certified. Its interior carries a ripple of amplitude
  0.0065, below the tolerances `-g 0.007 -e 0.007`, that leaves cells close to
  folding: the least-squares fit of a coarsened region folds them (LS alone ends at
  1994 functions, 787 of its 789 rejections for regularity), the linear stage
  removes every such fold (349 functions, all 120 rejections for the tolerances),
  and the non-linear stage is never invoked. The earlier input `triangle_wavy_L5.xml`
  (two folded cells in p0) and its runs `triangle_*` (3169, 581, 581 at `-g 0.01
  -e 0.01`) are kept for reference.
- The rotor runs of Table 5 use `-g 0.03 -e 0.0002`; with `-e 0.03` (`rotor_global`,
  `rotor_lambda0`, `rotor_lambda2`) every coarsening is accepted and all three settings end
  in the coarsest space. The cells and the fitting domain of the first candidate (Figure 8)
  are written by `--dump_first <file>`.

The times in the logs are those of single sequential runs on a Windows 11
laptop on mains power; they depend on the machine and its power state. The
final numbers of MPBES functions do not.

## Program options

`unrefinement_example --help` lists all options. The ones used above:

| Option | Meaning | Default |
|---|---|---|
| `-f` | input XML (multipatch of THB-spline patches with topology) | |
| `-g`, `-e` | tolerances epsilon_global (C2) and epsilon_feature (C3) | 1e6 (unconstrained) |
| `-F` | feature set F as `patch:side` items (side `u0`/`u1`/`v0`/`v1`, west/east/south/north or 1..4) | the whole outer boundary |
| `--interfaces` | measure (C3) on the interfaces as well | off |
| `-s` | cell selection `grendas` or `lex` | grendas |
| `-l` | locality parameter lambda of the fitting subdomain; -1 is global fitting | -1 |
| `-x` | factor f of the q_p = max(5, ceil(f sqrt(n_p))) points per direction of the grid of (C2) | 1 |
| `-p` | sample points of the fit: `gauss` (Gauss points of the elements) or `grid` (the q_p x q_p grid) | gauss |
| `-t` | stages: 1 = LS, 2 = LS+LO, 3 = LS+LO+NLO | 3 |
| `-o`, `--noplot` | output prefix; no Paraview output | |

Further options, not used in the runs (their defaults are the values of the paper):

| Option | Meaning | Default |
|---|---|---|
| `--b_linear`, `--b_nonlinear` | constant b of the adaptive weights of the linear and the non-linear stage | 2e-4, 2e-3 |
| `--beta_cap` | cap of n_-/10 in the factor beta = 1.35 (1 + min(cap, n_-/10)) of the weights | 10 |
| `--weight_floor` | floor of every weight as a fraction of beta b | 0.25 |
| `--scale_uniformity`, `--scale_length`, `--scale_orthogonality`, `--scale_skewness`, `--scale_eccentricity` | non-linear stage: factor on the weight of the functional after the weight rule (0 switches it off) | 1 |
| `--dump_first <file>` | writes the cells and the fitting domain of the first candidate (Figure 8) | |

## License

Mozilla Public License 2.0, as G+Smo (see `LICENSE`).
