"""
    WannierBSEKit

Umbrella package for Wannier-based Bethe--Salpeter equation (BSE) calculations.

`WannierBSEKit` contains no physics implementation of its own. It re-exports two
packages at exactly pinned, mutually compatible versions:

- [`WannierInterAction`](@ref): Coulomb interactions between Wannier pairs and
  mirror corrections of their long-range part.
- [`WannierExcitonModel`](@ref): tight-binding models, BSE Hamiltonians, exciton
  band structures and optical spectra.

`using WannierBSEKit` therefore provides the complete public API of both packages.
The version pins are part of the design: one release of `WannierBSEKit` identifies
one tested combination of the underlying packages.

# Getting help

The public API is documented with docstrings. In the Julia REPL press `?` to enter help
mode and type a name, for example

    help?> BSE
    help?> int
    help?> MirrorCorrection

or write `?BSE` directly at the `julia>` prompt. The package `README.md` contains the
installation instructions, worked examples and a comprehensive test run.

See also [`WannierInterAction`](@ref) and [`WannierExcitonModel`](@ref).
"""
module WannierBSEKit

using Reexport

@reexport using WannierInterAction
@reexport using WannierExcitonModel

end
