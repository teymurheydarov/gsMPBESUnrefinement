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
| `unittests/` | unit tests of the module (20 tests) |
| `filedata/` | `two_patches_refined.xml`, the input of the unit tests (G+Smo adds this folder to its file search path) |
| `data/geometries/` | the four input parameterizations of the paper |
| `runs/<run>/` | the logs of the 22 runs of Section 7 (see below) |
| `run_all.ps1` | reruns all experiments and compares the results with `runs/` (Windows) |

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
`-Only run1,run2` selects runs. The whole set takes a few minutes.

On other systems, every run can be started with the command line in the first
line of `runs/<run>/cmd_output.txt` (paths relative to the module folder).

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

| Run | Input | Paper | MPBES functions |
|---|---|---|---|
| `yeti_eg001_ef001`, `yeti_eg003_ef003`, `yeti_eg006_ef006` | Yeti | Table 2, row 1; Figure 2 | 239, 157, 122 |
| `yeti_nof_ef010`, `yeti_nof_ef030`, `yeti_nof_ef060` | Yeti | Table 2, row 2 | 239, 157, 122 |
| `yeti_F_ef010`, `yeti_F_ef030`, `yeti_F_ef060` | Yeti | Table 2, row 3; Figure 3 | 171, 128, 117 |
| `rotor_eg0001_ef0001`, `rotor_nog_ef0001` | rotor | Section 7.1; Figure 4 | 466, 251 |
| `mask_folded_ls_only`, `mask_folded_lo_only`, `mask_folded_nlo_full` | mask | Table 3; Figure 5 | 1089, 208, 130 |
| `joystick_grendas`, `joystick_lex` | joystick | Table 4; Figure 6 | 464, 472 |
| `joystick_global`, `joystick_lambda0`, `joystick_lambda1` | joystick | Section 7.2.3 | 464, 464, 464 |
| `rotor_global`, `rotor_lambda0`, `rotor_lambda2` | rotor | Table 5; Figures 7, 8 | 168, 168, 168 |

Notes on the settings:

- An unconstrained tolerance is the program's default `1e6`. The runs
  `yeti_nof_*` and `yeti_F_*` use `-g 1.0`, which never binds on the Yeti
  (the largest global error of a regular candidate in these runs is 0.086).
- The mask input is itself not regular (the log reports 6 cells not certified
  for the initial parameterization), so with LS alone no coarsening is
  accepted.
- The run time of `mask_folded_lo_only` exceeds that of
  `mask_folded_nlo_full`: without the non-linear stage most candidates are
  rejected (112 of 216, against 16 of 51), and every rejection leads to a
  further candidate.

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
| `-x` | factor f of the q_p = max(5, ceil(f sqrt(n_p))) sample points per direction | 1 |
| `-t` | stages: 1 = LS, 2 = LS+LO, 3 = LS+LO+NLO | 3 |
| `-o`, `--noplot` | output prefix; no Paraview output | |

## License

Mozilla Public License 2.0, as G+Smo (see `LICENSE`).
