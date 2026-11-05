using Frost

length(ARGS) <= 1 ||
    throw(ArgumentError("usage: julia ../run.jl [N]"))
N = isempty(ARGS) ? 16 : parse(Int, only(ARGS))

J1 = 1.0
J2 = 0.55

include(joinpath(@__DIR__, "run", "generate_lattice.jl"))
generate_j1j2_chain(
    N;
    J1,
    J2,
    output_dir = joinpath(@__DIR__, "run"),
)

results_file = "results.dat"
Frost.xtrg(N = N, D = 200, output_path = results_file)

include(joinpath(@__DIR__, "run", "plot_results.jl"))
plot_thermodynamics(results_file)
