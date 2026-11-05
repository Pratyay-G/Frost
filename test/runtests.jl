using DelimitedFiles
using LinearAlgebra
using Test
using Frost

include(joinpath(@__DIR__, "..", "run", "generate_lattice.jl"))
include(joinpath(@__DIR__, "..", "run", "plot_results.jl"))

@testset "Frost package" begin
    @test Frost.tantrg === Frost.TanTRG.main
    @test Frost.xtrg === Frost.XTRG.xtrg_main
    @test !isdefined(Frost.TanTRG, :xtrg_main)
    @test !isdefined(Frost.XTRG, :main)
    @test isdefined(Frost, :tantrg)
    @test isdefined(Frost, :xtrg)
    @test isfile(Frost.DEFAULT_LATTICE_INFO_PATH)
    @test_throws ArgumentError Frost.tantrg(N = 0)
end

@testset "XTRG thermal results" begin
    mktempdir() do output_dir
        lattice_path = joinpath(output_dir, "lattice_info.dat")
        results_path = joinpath(output_dir, "results.dat")
        open(lattice_path, "w") do io
            println(io, "site1\tsite2\tJ12")
            println(io, "1\t2\t1.0")
            println(io, "2\t3\t1.0")
        end

        Frost.xtrg(
            N = 3,
            D = 4,
            Nstep = 1,
            β0 = 1e-4,
            nz = 1,
            nmax = 1,
            lattice_info_path = lattice_path,
            output_path = results_path,
        )

        results = readdlm(results_path, Float64)
        @test size(results) == (1, 3)
        @test results[1, 1] == 2e-4
        @test all(isfinite, results)
        plot_path = joinpath(output_dir, "xtrg_properties.png")
        @test length(plot_thermodynamics(results_path, plot_path).plot.subplots) == 2
        @test isfile(plot_path)
    end
end

@testset "Thermal property plots" begin
    mktempdir() do output_dir
        results_path = joinpath(output_dir, "results.dat")
        plot_path = joinpath(output_dir, "thermal_properties.png")

        writedlm(results_path, [0.5 -1.0 4.0; 2.0 -0.5 1.0])
        result = plot_thermodynamics(results_path, plot_path)
        @test result.temperature == [0.5, 2.0]
        @test length(result.plot.subplots) == 2
        @test result.specific_heat === nothing
        @test propertynames(result) == (:plot, :temperature, :specific_heat)
        @test isfile(plot_path)

        temperatures = [0.5, 1.0, 2.0, 4.0, 7.0]
        energy = temperatures .^ 3 .- 2 .* temperatures .^ 2 .+ 3 .* temperatures .+ 4
        writedlm(results_path, hcat(1 ./ temperatures, zeros(5), energy))
        result = plot_thermodynamics(results_path, plot_path)
        expected_specific_heat = 3 .* temperatures .^ 2 .- 4 .* temperatures .+ 3
        @test length(result.plot.subplots) == 3
        @test result.temperature == temperatures
        @test result.specific_heat ≈ expected_specific_heat

        writedlm(results_path, [0.5 -1.0 4.0])
        result = plot_thermodynamics(results_path, plot_path)
        @test result.temperature == [2.0]
        @test result.specific_heat === nothing

        writedlm(results_path, [0.5 -1.0; 2.0 -0.5])
        @test_throws ArgumentError plot_thermodynamics(results_path, plot_path)

        writedlm(results_path, [0.0 -1.0 0.1])
        @test_throws ArgumentError plot_thermodynamics(results_path, plot_path)

        writedlm(results_path, [1.0 -1.0 0.0; 1.0 -0.5 0.2])
        @test_throws ArgumentError plot_thermodynamics(results_path, plot_path)
    end
end

@testset "Small-chain thermal results against exact diagonalization" begin
    mktempdir() do output_dir
        lattice_path = joinpath(output_dir, "lattice_info.dat")
        results_path = joinpath(output_dir, "results.dat")
        open(lattice_path, "w") do io
            println(io, "site1\tsite2\tJ12")
            for site in 1:3
                println(io, site, '\t', site + 1, '\t', 1.0)
            end
            for site in 1:2
                println(io, site, '\t', site + 2, '\t', 0.55)
            end
        end
        writedlm(results_path, [99.0 99.0 99.0])

        Frost.tantrg(
            N = 4,
            D = 64,
            βmax = 4.0,
            lattice_info_path = lattice_path,
            output_path = results_path,
        )
        results = readdlm(results_path, Float64)

        hamiltonian = zeros(16, 16)
        bonds = [(site, site + 1, 1.0) for site in 1:3]
        append!(bonds, [(site, site + 2, 0.55) for site in 1:2])
        for state in 0:15, (site1, site2, coupling) in bonds
            spin1 = ((state >> (site1 - 1)) & 1) == 1 ? 0.5 : -0.5
            spin2 = ((state >> (site2 - 1)) & 1) == 1 ? 0.5 : -0.5
            hamiltonian[state + 1, state + 1] += coupling * spin1 * spin2
            if spin1 != spin2
                flipped = state ⊻ (1 << (site1 - 1)) ⊻ (1 << (site2 - 1))
                hamiltonian[flipped + 1, state + 1] += coupling / 2
            end
        end
        exact_energies = eigvals(Symmetric(hamiltonian))
        β = results[1, 1]
        boltzmann_weights = exp.(-β .* exact_energies)
        exact_U = sum(exact_energies .* boltzmann_weights) / sum(boltzmann_weights)
        exact_F = -log(sum(boltzmann_weights)) / β

        @test size(results, 2) == 3
        @test size(results, 1) > 2
        @test isapprox(results[1, 2], exact_F; atol = 1e-3)
        @test isapprox(results[1, 3], exact_U; atol = 5e-4)

    end
end

@testset "Open J1-J2 spin-chain generation" begin
    mktempdir() do output_dir
        lattice = generate_j1j2_chain(8; J1 = 1.5, J2 = 0.25, output_dir)
        info, header = readdlm(joinpath(output_dir, "lattice_info.dat"); header = true)

        @test vec(header) == ["site1", "site2", "J12"]
        @test size(info) == (13, 3)
        @test count(==(1.5), info[:, 3]) == 7
        @test count(==(0.25), info[:, 3]) == 6
        @test all(isinteger, info[:, 1:2])
        @test all(1 .<= info[:, 1:2] .<= 8)
        @test isfile(joinpath(output_dir, "lattice.png"))
        @test length(lattice.j1_bonds) == 7
        @test length(lattice.j2_bonds) == 6
        @test (1, 2) in lattice.j1_bonds
        @test (1, 3) in lattice.j2_bonds
        @test !((1, 8) in lattice.j1_bonds)
        @test !((1, 7) in lattice.j2_bonds)
        @test_throws ArgumentError generate_j1j2_chain(0; output_dir)
    end
end
