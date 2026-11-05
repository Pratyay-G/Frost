# Frost

<p align="center">
  <img src="logo.png" alt="Frost logo" width="600">
</p>

Frost calculates thermal properties of finite quantum spin systems with two
independent algorithms: TanTRG and XTRG. They are available as separate
modules (`Frost.TanTRG` and `Frost.XTRG`) and through the package-level
`Frost.tantrg` and `Frost.xtrg` entry points.

## Setup

Instantiate the Julia environment from the repository root:

```julia
using Pkg
Pkg.instantiate()
```

## Generate a lattice and run the simulation

This package simulates the open-boundary one-dimensional J1-J2 spin chain.
Nearest-neighbor bonds have coupling J1 and next-nearest-neighbor bonds have
coupling J2. Run the simulation from the `run/` directory. It defaults to 28
sites, writes the bonds to `lattice_info.dat` (columns `site1 site2 J12`),
saves a site-labeled plot as `lattice.png`, and then runs the calculation.
Set `J1` and `J2` near the top of `../run.jl` to choose the couplings.

```sh
julia --project=.. ../run.jl
```

Pass a site count to simulate a chain of a different length, for example:

```sh
julia --project=.. ../run.jl 32
```

The run script writes β, free energy F, and internal energy U to
`energy_vs_beta_with_U.dat`, replacing results from the previous run. Run
`plot_results.jl` separately from the `run/`
directory at any time to create `thermal_properties.png`, which plots F and U
versus temperature T = 1/β, and specific heat C = dU/dT when at least four
distinct temperature points are available. Specific heat is estimated with a
four-point finite difference:

```sh
julia --project=.. plot_results.jl energy_vs_beta_with_U.dat thermal_properties.png
```

The plotting script requires a U column and throws an error if it is missing.
`run.jl` runs it automatically after the calculation.
The default MPS truncation cutoff is `1e-8` to reduce accumulated errors in
internal energy.
Other simulation parameters can be overridden by passing keywords to
`Frost.tantrg`, for example:

```julia
Frost.tantrg(N = 16, D = 200, βmax = 40.0, output_path = "results.dat")
```

The XTRG algorithm uses the same lattice input and writes β, F, and U
columns. Run it explicitly with:

```julia
Frost.xtrg(N = 16, D = 200, output_path = "results.dat")
```

Its sweep count, bond dimension, initial β, number of logarithmic starting
points, and number of doubling steps can be configured with `Nstep`, `D`,
`β0`, `nz`, and `nmax`, respectively.
