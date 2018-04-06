# ------------------------------------------------------------------
# Copyright (c) 2018, Júlio Hoffimann Mendes <juliohm@stanford.edu>
# Licensed under the ISC License. See LICENCE in the project root.
# ------------------------------------------------------------------

__precompile__()

module DirectGaussianSimulation

importall GeoStatsBase
using GeoStatsDevTools

using Variography

"""
    DirectGaussSim(var₁=>param₁, var₂=>param₂, ...)

Direct Gaussian simulation (a.k.a. LU simulation).

## Parameters

* `variogram` - theoretical variogram (default to GaussianVariogram())

### References

Alabert 1987. *The practice of fast conditional simulations through the
LU decomposition of the covariance matrix.*
"""
@simsolver DirectGaussSim begin
  @param variogram = GaussianVariogram()
end

function solve(problem::SimulationProblem, solver::DirectGaussSim)
  # sanity checks
  @assert keys(solver.params) ⊆ keys(variables(problem)) "invalid variable names in solver parameters"

  # retrieve problem domain
  pdomain = domain(problem)

  realizations = []
  for (var,V) in variables(problem)
    # get user parameters
    if var ∈ keys(solver.params)
      varparams = solver.params[var]
    else
      varparams = DirectGaussSimParam()
    end

    # determine variogram model
    γ = varparams.variogram

    # check stationarity
    @assert isstationary(γ) "variogram model must be stationary"

    # build covariance matrix
    C = γ.sill - pairwise(γ, pdomain, 1:npoints(pdomain))

    # Cholesky factorization
    L = chol(C)'

    # retrieve data locations in domain
    datalocs = [loc for (loc, dataloc) in datamap(problem, var)]

    if nworkers() > 1
      # generate realizations in parallel
      λ = _ -> solve_single(problem, var, L)
      varreals = pmap(λ, 1:nreals(problem))
    else
      # fallback to serial execution
      varreals = [solve_single(problem, var, L) for i=1:nreals(problem)]
    end

    push!(realizations, var => varreals)
  end

  SimulationSolution(domain(problem), Dict(realizations))
end

function solve_single(problem::AbstractProblem, var::Symbol, L::LowerTriangular)
  w = randn(size(L, 1))
  L*w
end

export DirectGaussSim

end
