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

  # retrieve problem info
  pdata = data(problem)
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

    # retrieve data locations in domain and data values
    datalocs = Vector{Int}()
    z₁ = Vector{V}()
    for (loc, dataloc) in datamap(problem, var)
      push!(datalocs, loc)
      push!(z₁, value(pdata, dataloc, var))
    end

    # retrieve simulation locations
    simulated = falses(npoints(pdomain))
    simulated[datalocs] = true
    simlocs = find(.!simulated)

    # covariance between simulation locations
    C₂₂ = sill(γ) - pairwise(γ, pdomain, simlocs)

    if isempty(datalocs)
      d₂  = zero(V)
      L₂₂ = chol(Symmetric(C₂₂))'
    else
      # covariance beween data locations
      C₁₁ = sill(γ) - pairwise(γ, pdomain, datalocs)
      C₁₂ = sill(γ) - pairwise(γ, pdomain, datalocs, simlocs)

      L₁₁ = chol(Symmetric(C₁₁))'
      B₁₂ = L₁₁ \ C₁₂
      A₂₁ = B₁₂'

      d₂ = A₂₁*(L₁₁ \ z₁)
      L₂₂ = chol(Symmetric(C₂₂ - A₂₁*B₁₂))'
    end

    if nworkers() > 1
      # generate realizations in parallel
      λ = _ -> solve_single(problem, var, z₁, d₂, L₂₂, datalocs, simlocs)
      varreals = pmap(λ, 1:nreals(problem))
    else
      # fallback to serial execution
      varreals = [solve_single(problem, var, z₁, d₂, L₂₂, datalocs, simlocs) for i=1:nreals(problem)]
    end

    push!(realizations, var => varreals)
  end

  SimulationSolution(domain(problem), Dict(realizations))
end

function solve_single(problem::AbstractProblem, var::Symbol,
                      z₁::Vector, d₂::Union{Real,Vector}, L₂₂::LowerTriangular,
                      datalocs::Vector{Int}, simlocs::Vector{Int})
  # retrieve problem info
  pdomain = domain(problem)
  V = variables(problem)[var]

  # allocate memory for result
  realization = Vector{V}(npoints(pdomain))

  # set hard data
  realization[datalocs] = z₁

  # simulate the rest
  w₂ = randn(size(L₂₂, 2))
  y₂ = d₂ + L₂₂*w₂
  realization[simlocs] = y₂

  realization
end

export DirectGaussSim

end
