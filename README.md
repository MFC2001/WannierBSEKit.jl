# WannierBSEKit.jl

[![Build Status](https://github.com/MFC2001/WannierBSEKit.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/MFC2001/WannierBSEKit.jl/actions/workflows/CI.yml?query=branch%3Amaster)

**An umbrella package for Wannier-based Bethe–Salpeter equation calculations.**

`WannierBSEKit.jl` is the entry point to a versioned software stack that computes excitonic states directly 
in an electronic Wannier basis. 
A single `using WannierBSEKit` exposes the complete public API of two scientific packages, pinned to mutually compatible versions:

- **[WannierInterAction.jl](https://github.com/MFC2001/WannierInterAction.jl)** — Coulomb interactions between
  Wannier pairs (`Rhok`, `int`) and mirror correction of their long-range part (`MirrorCorrection`);
- **[WannierExcitonModel.jl](https://github.com/MFC2001/WannierExcitonModel.jl)** — tight-binding models from
  Wannier90 output (`TightBindModel`), BSE Hamiltonians (`BSE`), exciton bands (`BAND`) and optical spectra
  (`absorption`, `optcond`).

The package also ships three worked examples — graphone, LiF and monolayer MoS2 — with their input files,
driver scripts and reference output data.

> **Scope.** This is an *umbrella* package: it contains no physics implementation of its own. It fixes one tested
> combination of the underlying packages and distributes the examples. It is **not registered in the General
> registry**; install it as described in §3.

## 1. Capabilities

- Direct, exchange and pair Coulomb interactions between Wannier functions for **2D slabs** and **3D bulk**, using
  bare or screened (Yambo, BerkeleyGW) dielectric matrices.
- Mirror correction of the periodic-image tail (`:gauss`, `:gauss_optim`, `:keldysh_gauss`, `:keldysh_gauss_optim`).
- Spin-unpolarized, spin-polarized and noncollinear (spinor, spin–orbit) BSE Hamiltonians, with an optional scissor
  shift.
- Exciton band structures with `LinearAlgebra` or `KrylovKit`; optical conductivity (`optcond`), absorption
  (`absorption`); exciton wavefunctions and exciton Wannier functions.
- Readers and writers for Quantum ESPRESSO, Wannier90, WannierBerri, BerkeleyGW and Yambo files.

## 2. Documentation and the built-in help system

The reference documentation is the set of docstrings in the installed source, available offline through **help
mode**. At the `julia>` prompt press `?`; the prompt becomes `help?>`, where a name can be typed:

```julia
julia> using WannierBSEKit

julia> ?
help?> BSE
```

Backspace on an empty `help?>` line (or `Ctrl+C`) returns to the normal prompt. 
The top-level entries are `?WannierBSEKit` (overview of the bundle and this instruction),
`?WannierInterAction` and `?WannierExcitonModel`.

| Task | Help |
|---|---|
| Real space | `?Cell`, `?Lattice`, `?ReducedCoordinates`, `?CartesianCoordinates` |
| Reciprocal space | `?MonkhorstPack`, `?RedKgrid`, `?Kline` |
| Wavefunctions and screening | `?Rhok`, `?QE_wannier`, `?epsmat_bgw`, `?epsmat_yambo` |
| Interactions and correction | `?int`, `?MirrorCorrection`, `?PairInterAction` |
| Tight-binding model and BSE model | `?TightBindModel`, `?BSE` |
| Exciton bands and optical spectra | `?BAND`, `?OpticalTransitionData`, `?absorption`, `?optcond` |

## 3. Installation

`WannierBSEKit` requires Julia 1.10 or newer. Its dependencies are all registered in the General registry, so only
`WannierBSEKit` itself has to be installed explicitly; Julia installs the dependencies automatically at the
versions specified in `Project.toml`.

### 3.1 On a networked machine

Install the package from its repository with `Pkg.add`:

```julia
using Pkg
Pkg.add(url = "https://github.com/MFC2001/WannierBSEKit.jl")
```

Alternatively, clone the repository and install the local copy with `Pkg.develop`:

```julia
Pkg.develop(path = "/path/to/WannierBSEKit")
```

The example scripts additionally use `JLD2` and `KrylovKit`; install them into the same environment:

```julia
Pkg.add(["JLD2", "KrylovKit"])
```

### 3.2 On an offline machine

1. On a networked machine with the same operating system and the same Julia version, install `WannierBSEKit` as
   described in §3.1 — preferably with `Pkg.add`, so that the package source is kept inside `~/.julia` — together
   with the other Julia packages required by the workflow (for example `KrylovKit` and `JLD2`). This populates the
   `~/.julia` directory.
2. Copy `~/.julia` from the networked machine and use it to replace `~/.julia` on the offline machine.

The complete depot is transferred, so the offline machine needs no package-manager access. The operating system and
Julia version must match the networked machine, because the depot contains platform-specific precompiled caches.

### 3.3 External software used by the examples

The Julia package itself requires no external program. The worked examples additionally invoke Quantum ESPRESSO,
Wannier90, BerkeleyGW, Yambo and WannierBerri; the corresponding commands and input files are given in the
`examples/` directories.

## 4. Repository and file layout

```text
WannierBSEKit/
├── Project.toml                 package metadata, exact version pins
├── LICENSE                      MIT
├── src/WannierBSEKit.jl         module: docstring and re-exports only
├── test/
├── .github/
└── examples/
    ├── Graphone/                2D ferromagnetic graphone (Yambo screening)
    ├── LiF/                     3D bulk LiF (BerkeleyGW screening)
    └── MoS2/                    2D monolayer MoS2 (noncollinear, BerkeleyGW)
```

Apart from the module file, all functionality is in the pinned dependencies; there are no other source files.

## 5. Running the examples

### 5.1 How the `script` files are used

Each material directory contains `script`, a **command checklist** rather than a self-contained shell program. It
lists the commands in execution order, grouped by the step directory in the comment (`# 01-scf`, ...). Each program
block is executed **inside its step directory**; the `cp`/`ln`/`mv` preparation lines are executed from the material
directory. Before running, set the prefixes used by the scripts:

```bash
export QEbin=/path/to/qe/bin
export bgwbin=/path/to/berkeleygw/bin
export yambobin=/path/to/yambo/bin
export juliabin=/path/to/julia/bin
```

and replace `pseudo_dir` in every `pw.in`.

### 5.2 Practical notes

- **Use Graphone for a quick check.** If you only want to verify that the whole workflow runs on your machine,
  start with [`examples/Graphone`](examples/Graphone): it is the smallest and fastest of the three cases.
- **The numerical parameters are production values.** Apart from file paths (and the executable prefixes above),
  every parameter in the template inputs — cutoffs, k-point meshes, band windows, scissor shifts, broadening,
  interaction cutoffs, etc. — is the production value used to obtain the reference results in the example folders.
  They are not arbitrary placeholders: the final exciton bands and absorption or conductivity spectra are reproduced
  only with these settings, and changing them changes the results.
- **Split the screening step of LiF and MoS2.** The dielectric-function step has many k-points, and a single
  `epsilon.cplx.x` job that treats all of them may fail. Run several jobs, each covering a subset of the k-points
  (the `begin qpoints` block in `epsilon.inp` can be restricted accordingly), then merge the resulting dielectric
  functions with BerkeleyGW's `$bgwbin/epsmat_hdf5_merge.py`. The multiple output files in the epsilon step come
  from these partial jobs.

### 5.3 What each case demonstrates

- **Graphone** — 2D ferromagnet with two spin channels, Yambo screening, spin-resolved BSE band
  (`07-bse/bseband.dat`) and one-particle bands (`eleband/`).
- **LiF** — 3D bulk with BerkeleyGW screening, WannierBerri Hamiltonian and a GW scissor shift.
- **MoS2** — noncollinear monolayer with spin–orbit coupling, 22 spinor orbitals, `channel = 1`.

The unnumbered folders (`eleband`, `excwave`, `excband_Gamma`, `yambo-bse`, `yambo-gw`) are independent
post-processing or benchmark tasks; they use the output of the main workflow but do not modify it.
