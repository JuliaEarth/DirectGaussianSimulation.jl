# ------------------------------------------------------------------
# Copyright (c) 2018, Júlio Hoffimann Mendes <juliohm@stanford.edu>
# Licensed under the ISC License. See LICENCE in the project root.
# ------------------------------------------------------------------

__precompile__()

module DirectGaussianSimulation

importall GeoStatsBase
using GeoStatsDevTools

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
end

export DirectGaussSim

end
