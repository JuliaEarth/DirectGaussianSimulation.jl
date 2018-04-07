using GeoStats
using DirectGaussianSimulation
using Plots; gr()
using Base.Test
using VisualRegressionTests

# setup GR backend for Travis CI
ENV["GKSwstype"] = "100"
ENV["PLOTS_TEST"] = "true"

# list of maintainers
maintainers = ["juliohm"]

# environment settings
ismaintainer = "USER" ∈ keys(ENV) && ENV["USER"] ∈ maintainers
istravislinux = "TRAVIS" ∈ keys(ENV) && ENV["TRAVIS_OS_NAME"] == "linux"
datadir = joinpath(@__DIR__,"data")

@testset "DirectGaussianSimulation.jl" begin
  geodata = GeoDataFrame(DataFrames.DataFrame(x=[0.,25.,50.,75.,100.], variable=[0.,1.,0.,1.,0.]), [:x])
  domain = RegularGrid{Float64}(100)
  @testset "Conditional simulation" begin
    problem = SimulationProblem(geodata, domain, :variable, 3)

    srand(2018)
    solver = DirectGaussSim(:variable => @NT(variogram=SphericalVariogram(range=10.)))

    solution = solve(problem, solver)

    if ismaintainer || istravislinux
      function plot_cond_solution(fname)
        plot(solution, size=(1000,400))
        png(fname)
      end
      refimg = joinpath(datadir,"CondSimSol.png")

      @test test_images(VisualTest(plot_cond_solution, refimg), popup=!istravislinux) |> success
    end
  end

  @testset "Unconditional simulation" begin
    problem = SimulationProblem(domain, :variable => Float64, 3)

    srand(2018)
    solver = DirectGaussSim(:variable => @NT(variogram=SphericalVariogram(range=10.)))

    solution = solve(problem, solver)

    if ismaintainer || istravislinux
      function plot_uncond_solution(fname)
        plot(solution, size=(1000,400))
        png(fname)
      end
      refimg = joinpath(datadir,"UncondSimSol.png")

      @test test_images(VisualTest(plot_uncond_solution, refimg), popup=!istravislinux) |> success
    end
  end
end
