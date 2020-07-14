using DirectGaussianSimulation
using GeoStatsBase
using Variography
using Plots, VisualRegressionTests
using Test, Pkg, Random

# workaround GR warnings
ENV["GKSwstype"] = "100"

# environment settings
islinux = Sys.islinux()
istravis = "TRAVIS" ∈ keys(ENV)
datadir = joinpath(@__DIR__,"data")
visualtests = !istravis || (istravis && islinux)
if !istravis
  Pkg.add("Gtk")
  using Gtk
end

@testset "DirectGaussianSimulation.jl" begin
  @testset "Conditional simulation" begin
    sdata = georef((z=[0.,1.,0.,1.,0.],), [0. 25. 50. 75. 100.])
    sdomain = RegularGrid(100)

    problem = SimulationProblem(sdata, sdomain, :z, 2)
    solver = DirectGaussSim(:z => (variogram=SphericalVariogram(range=10.),))

    Random.seed!(2018)
    solution = solve(problem, solver)

    if visualtests
      @plottest plot(solution) joinpath(datadir,"CondSimSol.png") !istravis
    end
  end

  @testset "Unconditional simulation" begin
    sdomain = RegularGrid(100)

    problem = SimulationProblem(sdomain, :z => Float64, 2)
    solver = DirectGaussSim(:z => (variogram=SphericalVariogram(range=10.),))

    Random.seed!(2018)
    solution = solve(problem, solver)

    if visualtests
      @plottest plot(solution) joinpath(datadir,"UncondSimSol.png") !istravis
    end
  end

  @testset "Cosimulation" begin
    sdata = georef((z=[0.,1.,0.,1.,0.], y=[0.,1.,0.,1.,0.]), [0. 25. 50. 75. 100.])
    sdomain = RegularGrid(500)

    problem = SimulationProblem(sdomain, (:z=>Float64,:y=>Float64), 1)
    solver = DirectGaussSim(:z => (variogram=SphericalVariogram(range=10.),),
                            :y => (variogram=GaussianVariogram(range=10.),),
                            (:z,:y) => (correlation=0.95,))

    Random.seed!(2020)
    solution = solve(problem, solver)

    if visualtests
      @plottest plot(solution) joinpath(datadir,"CoSimSol.png") !istravis
    end
  end
end
