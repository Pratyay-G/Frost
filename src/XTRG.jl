module XTRG

using DelimitedFiles
using ITensorMPS
using ITensors
using LinearAlgebra

export xtrg_main

include("sweeps.jl")

const DEFAULT_LATTICE_INFO_PATH = normpath(joinpath(@__DIR__, "..", "run", "lattice_info.dat"))

function xtrg_main(;
    N::Union{Nothing,Integer} = nothing,
    D::Integer = 20,
    Nstep::Integer = 4,
    β0::Real = 1e-6,
    nz::Integer = 4,
    nmax::Integer = 28,
    lattice_info_path::AbstractString = DEFAULT_LATTICE_INFO_PATH,
    output_path::AbstractString = "results.dat",
)
    isnothing(N) || N > 0 || throw(ArgumentError("N must be positive"))
    D > 0 || throw(ArgumentError("D must be positive"))
    Nstep > 0 || throw(ArgumentError("Nstep must be positive"))
    isfinite(β0) && β0 > 0 || throw(ArgumentError("β0 must be finite and positive"))
    nz > 0 || throw(ArgumentError("nz must be positive"))
    nmax > 0 || throw(ArgumentError("nmax must be positive"))

    BLAS.set_num_threads(1)
    ITensors.Strided.set_num_threads(1)
    ITensors.enable_threaded_blocksparse()

    bonds, header = readdlm(lattice_info_path; header = true)
    vec(header) == ["site1", "site2", "J12"] ||
        throw(ArgumentError("lattice info header must be: site1 site2 J12"))
    size(bonds, 2) == 3 ||
        throw(ArgumentError("lattice info must contain three columns: site1 site2 J12"))
    all(isinteger, bonds[:, 1:2]) ||
        throw(ArgumentError("lattice info site indices must be integers"))

    inferred_N = isempty(bonds) ? nothing : Int(maximum(bonds[:, 1:2]))
    n_sites = isnothing(N) ? inferred_N : N
    isnothing(n_sites) &&
        throw(ArgumentError("N must be provided when lattice info contains no bonds"))
    n_sites > 0 || throw(ArgumentError("N must be positive"))

    sites = siteinds("S=1/2", n_sites; conserve_qns = true)
    osH = OpSum()
    for (site1, site2, coupling) in eachrow(bonds)
        1 <= site1 <= n_sites && 1 <= site2 <= n_sites ||
            throw(ArgumentError("lattice info contains a site index outside 1:N"))
        isfinite(coupling) ||
            throw(ArgumentError("lattice info coupling values must be finite"))
        osH += (coupling, "Sz", Int(site1), "Sz", Int(site2))
        osH += (0.5 * coupling, "S+", Int(site1), "S-", Int(site2))
        osH += (0.5 * coupling, "S-", Int(site1), "S+", Int(site2))
    end

    ham = MPO(osH, sites)
    osI = OpSum()
    osI += ("Id", 1)
    identity_mpo = MPO(osI, sites)

    open(output_path, "w") do io
        for seed_index in 0:(nz - 1)
            β_current = β0 * 2.0^(seed_index / nz)
            rho = identity_mpo - β_current * ham
            log_partition_scale = 0.0

            for _ in 1:nmax
                β_current *= 2.0
                println("β = ", β_current)

                C = dag(prime(prime(rho)))
                orthogonalize!(C, n_sites)

                for _ in 1:Nstep
                    C, rho, log_partition_scale =
                        left2right(rho, C, n_sites, D, sites, log_partition_scale)
                    C, rho, log_partition_scale =
                        right2left(rho, C, n_sites, D, sites, log_partition_scale)
                end

                for site in 1:n_sites
                    rho[site] =
                        (dag(C[site]) * delta(dag(sites[site]'''), sites[site]')) *
                        delta(sites[site]'', dag(sites[site]))
                end

                partition_function = 1.0
                for site in 1:n_sites
                    partition_function *=
                        rho[site] * delta(dag(sites[site]'), sites[site])
                end
                Z = real(partition_function[])
                isfinite(Z) && Z > 0 ||
                    throw(ArgumentError("XTRG produced an invalid partition function"))

                energy_numerator = 1.0
                for site in 1:n_sites
                    energy_numerator = (energy_numerator * rho[site]) * dag(ham[site])
                end
                U = real(energy_numerator[] / Z)

                log_partition_scale *= 2.0
                F = -(log(Z) + log_partition_scale) / β_current
                writedlm(io, [β_current F U])
                flush(io)
            end
        end
    end

    println("Done")
    return nothing
end

end
