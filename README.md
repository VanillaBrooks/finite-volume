# Sod Shock Tube Solver


## Architecture

The repository contains both a julia and a rust implementation. 

`./src/SodsShockTube.jl` handles main solver loop functionality, including initializing variables
and matricies

`./src/state.jl` contains the definitions of the initial conditions and a `State` struct containing
information on density, pressure, velocity, entropy, etc.

`./src/io_utils.jl` has helper functions for writing output HDF5 files

`./src/solver_methods.jl` implements the various solver methods, generic to the argument to `main` in 
`SodsShockTube.jl`

`./src/vector_calculations.jl` has helper methods for calculating the `State` variables from `q`, calculating
`q` from the state variables, and calculating `F` from the state variables.


### Rust differences

The only difference between the rust and julia code is that `SodsShockTube.jl` is equivalent to `main.rs`

## Running

First, ensure you have a `./results` folder for the hdf5 files to be written to


To run the julia code:

```
julia ./src/SodsShockTube.jl
```

To run the rust code:

```
cargo r --release
```

## Plotting

Plotting is done through `Pluto` notebooks. If julia is installed:

```
julia
>]
(@v1.7) pkg> add Pluto
```

and then, from a julia REPL:

```
using Pluto
Pluto.run()
```

Since julia stores arrays in column major order, and the rust ndarray library stores
them in row major order, you must transpose the arrays read from HDF5 files. Set
`IS_RUST = true` for rust outputs, and `IS_RUST = false` for julia outputs.
