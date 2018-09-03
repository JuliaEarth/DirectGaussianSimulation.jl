using GeoStatsBase
using GeoStatsDevTools
using DirectGaussianSimulation
using Plots; gr()
using VisualRegressionTests
using Test, Pkg, Random

# list of maintainers
maintainers = ["juliohm"]

# environment settings
istravis = "TRAVIS" ∈ keys(ENV)
ismaintainer = "USER" ∈ keys(ENV) && ENV["USER"] ∈ maintainers
datadir = joinpath(@__DIR__,"data")

if ismaintainer
  Pkg.add("Gtk")
  using Gtk
end

@testset "DirectGaussianSimulation.jl" begin
  geodata = PointSetData(Dict(:z => [0.,1.,0.,1.,0.]), [0. 25. 50. 75. 100.])
  domain = RegularGrid{Float64}(100)

  @testset "Conditional simulation" begin
    problem = SimulationProblem(geodata, domain, :z, 2)

    Random.seed!(2018)
    solver = DirectGaussSim(:z => (variogram=SphericalVariogram(range=10.),))

    solution = solve(problem, solver)

    if ismaintainer || istravis
      function plot_cond_solution(fname)
        plot(solution, size=(1000,400))
        png(fname)
      end
      refimg = joinpath(datadir,"CondSimSol.png")

      @test test_images(VisualTest(plot_cond_solution, refimg), popup=!istravis, tol=0.1) |> success
    end
  end

  @testset "Unconditional simulation" begin
    problem = SimulationProblem(domain, :z => Float64, 2)

    Random.seed!(2018)
    solver = DirectGaussSim(:z => (variogram=SphericalVariogram(range=10.),))

    solution = solve(problem, solver)

    if ismaintainer || istravis
      function plot_uncond_solution(fname)
        plot(solution, size=(1000,400))
        png(fname)
      end
      refimg = joinpath(datadir,"UncondSimSol.png")

      @test test_images(VisualTest(plot_uncond_solution, refimg), popup=!istravis, tol=0.1) |> success
    end
  end
end
