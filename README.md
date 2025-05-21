# Rosenbluth

[![Build Status](https://github.com/tomrobertsof/Rosenbluth.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/tomrobertsof/Rosenbluth.jl/actions/workflows/CI.yml?query=branch%3Amain)

Rosenbluth.jl is a Julia library designed to provide a comprehensive suite of sampling algorithms within the Rosenbluth family. This library aims to offer efficient, easy-to-use implementations of these algorithms for applications in statistical physics, Monte Carlo simulations, and other fields requiring random sampling techniques.

## Features

- Implementation of the original Rosenbluth and Rosenbluth (RR) algorithm for chain molecule simulations.
- Advanced sampling techniques including the pruned-enriched Rosenbluth method (PERM) and the generalized atmospheric Rosenbluth method (GARM)
- Minimal, extensible design allowing for easy use on any model

## Installation

To install Rosenbluth.jl, you can use the Julia package manager. Open the Julia command line and run:

```julia
using Pkg
Pkg.add("https://github.com/tomrobertsOF/Rosenbluth.jl")
```

## Usage
The easiest way to get started is using the pre-defined models in the `Models` sub-module, for example:
```julia
using Rosenbluth, Rosenbluth.Models

# Runs 1000 tours of site trees up to size 50 with normal pruning/enrichment
sample(Models.SiteTree, 50, 1000; prune_enrich_method=PruneEnrichMethod.STANDARD )

```

In order to use your own model, define it as a type deriving from `RosenbluthSampleable` and implement methods for `atmosphere`, `size` and `grow!`.

```julia
using Rosenbluth
struct Model <: RosenbluthSampleable
    ...
end
function Rosenbluth.atmosphere(model::Model)::Int
    ...
end
function Rosenbluth.size(model::Model)::Int
    ...
end
function Rosenbluth.grow!(model::Model)
    ...
end
```

You can then call one of the sampling functions, for example:
```julia
sample(Model, max_size, num_runs)
```


## Contributing
Contributions to Rosenbluth.jl are very welcome! If you have a suggestion for an improvement or have found a bug, please open an issue on our GitHub repository. If you'd like to contribute code, please submit a pull request.

## License
Rosenbluth.jl is licensed under the MIT License. See the LICENSE file for more details.

## TODO

- [ ] Finalize the core API for the sampling algorithms.
- [x] Implement the basic Rosenbluth and Rosenbluth (RR) sampling algorithm.
- [x] Add support for the pruned-enriched Rosenbluth method (PERM).
- [ ] Create a suite of examples demonstrating the use of each algorithm.
- [x] Set up a testing framework and add initial unit tests.
- [ ] Optimize the performance of the algorithms for large-scale simulations.
- [x] Allow smart optimizations by flags or properties of the model
- [ ] Document the API and provide a getting started guide for new users.
- [x] Set up continuous integration (CI) to automate testing and ensure code quality.

