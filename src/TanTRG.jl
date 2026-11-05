module TanTRG

using DelimitedFiles
using ITensorMPS
using ITensors
using LinearAlgebra

export choi_mpo_to_mps, main

include("choi.jl")

const DEFAULT_LATTICE_INFO_PATH = normpath(joinpath(@__DIR__, "..", "run", "lattice_info.dat"))

function main(;
    calculate_U::Bool = true,
    N::Union{Nothing,Integer} = nothing,
    D::Integer = 400,
    Nstep::Integer = 1,
    cut_off::Real = 1e-8,
    β0::Real = 0.000001,
    βmax::Real = 20.0,
    lattice_info_path::AbstractString = DEFAULT_LATTICE_INFO_PATH,
    output_path::AbstractString = "energy_vs_beta.dat",
)
    isnothing(N) || N > 0 || throw(ArgumentError("N must be positive"))
    D > 0 || throw(ArgumentError("D must be positive"))
    Nstep > 0 || throw(ArgumentError("Nstep must be positive"))
    cut_off >= 0 || throw(ArgumentError("cut_off must be nonnegative"))
    β0 >= 0 || throw(ArgumentError("β0 must be nonnegative"))
    βmax > 2 * β0 || throw(ArgumentError("βmax must be greater than 2β0"))

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
    idTensor = MPO(osI, sites)
    rho_0 = idTensor - β0 * ham

    ϕ, mod_ham = choi_mpo_to_mps(rho_0, ham)

    β_current = 2.0 * β0
    numer = inner(ϕ', mod_ham, ϕ)
    Z = norm(ϕ)^2
    U = real(numer / Z)
    println("Starting at ", β_current, " and U = ", U)
    flush(stdout)

    logfact = log(norm(ϕ))
    ϕ = normalize(ϕ)

    write(output_path, "")

    while β_current < βmax
        remaining_β = βmax - β_current
        proposed_δβ = max(0.01, 0.01 * β_current)
        reaches_βmax = 2 * proposed_δβ >= remaining_β
        δβ = reaches_βmax ? remaining_β / 2 : proposed_δβ
        ϕ = tdvp(
            mod_ham,
            -δβ,
            ϕ;
            nsweeps = Nstep,
            maxdim = D,
            cutoff = cut_off,
            updater_kwargs = (; tol = 1e-6, krylovdim = 20),
            normalize = false,
            reverse_step = true,
            outputlevel = 1,
        )

        β_next = β_current + 2.0 * δβ
        if reaches_βmax || isapprox(β_next, βmax; atol = eps(βmax), rtol = 0)
            β_current = βmax
        else
            β_current = β_next
        end
        Z = norm(ϕ)^2
        F = (-log(Z) - 2.0 * logfact) / β_current
        println(β_current, " ", F)

        open(output_path, "a") do io
            if calculate_U
                numer = inner(ϕ', mod_ham, ϕ)
                U = real(numer / Z)
                writedlm(io, [β_current F U])
            else
                writedlm(io, [β_current F])
            end
        end
        flush(stdout)

        logfact += log(norm(ϕ))
        ϕ = normalize(ϕ)
    end

    println("Done")
    return nothing
end

end
