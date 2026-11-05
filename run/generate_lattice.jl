using Plots

function generate_j1j2_chain(
    N::Integer = 28;
    J1::Real = 1.0,
    J2::Real = 1.0,
    output_dir::AbstractString = @__DIR__,
)
    N > 0 || throw(ArgumentError("N must be positive"))

    j1_bonds = [(site, site + 1) for site in 1:(N - 1)]
    j2_bonds = [(site, site + 2) for site in 1:(N - 2)]

    mkpath(output_dir)
    open(joinpath(output_dir, "lattice_info.dat"), "w") do io
        println(io, "site1\tsite2\tJ12")
        for (bonds, coupling) in ((j1_bonds, J1), (j2_bonds, J2))
            for (site1, site2) in bonds
                println(io, site1, '\t', site2, '\t', coupling)
            end
        end
    end

    p = plot(
        aspect_ratio = :equal,
        xlims = (0.5, N + 0.5),
        ylims = (-0.75, 0.75),
        xlabel = "Site",
        ylabel = "",
        title = "Open J1-J2 spin chain (N = $N)",
        grid = false,
        yticks = false,
    )

    for (bonds, color, style, height) in (
        (j1_bonds, :black, :solid, -0.08),
        (j2_bonds, :red, :dash, 0.3),
    )
        for (i, j) in bonds
            t = range(0, 1; length = 30)
            xs = i .+ (j - i) .* t
            ys = 4 * height .* t .* (1 .- t)
            plot!(p, xs, ys; color, linewidth = 1.5, linestyle = style, label = nothing)
        end
    end

    scatter!(
        p,
        1:N,
        zeros(N);
        color = :white,
        markerstrokecolor = :black,
        markersize = 6,
        label = "Sites",
    )
    annotate!(p, [(site, -0.2, text(string(site), 7, :black)) for site in 1:N])
    plot!(p, [NaN], [NaN]; color = :black, linewidth = 1.5, label = "J1 bonds")
    plot!(p, [NaN], [NaN]; color = :red, linewidth = 1.5, linestyle = :dash, label = "J2 bonds")
    savefig(p, joinpath(output_dir, "lattice.png"))

    return (j1_bonds, j2_bonds, plot = p)
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) <= 1 ||
        throw(ArgumentError("usage: julia generate_lattice.jl [N]"))
    N = isempty(ARGS) ? 28 : parse(Int, only(ARGS))
    generate_j1j2_chain(N)
end
