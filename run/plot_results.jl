using DelimitedFiles
using Plots

function four_point_derivative(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    length(x) == length(y) || throw(ArgumentError("x and y must have the same length"))
    length(x) >= 4 || throw(ArgumentError("four-point differentiation requires at least four points"))

    derivative = Vector{Float64}(undef, length(x))
    for i in eachindex(x)
        first_index = clamp(i - 1, firstindex(x), lastindex(x) - 3)
        indices = first_index:(first_index + 3)
        offsets = Float64.(x[indices] .- x[i])
        scale = maximum(abs, offsets)
        nodes = offsets ./ scale
        weights = zeros(4)

        for j in 1:4
            denominator = 1.0
            for k in 1:4
                k == j || (denominator *= nodes[j] - nodes[k])
            end

            numerator = 0.0
            for m in 1:4
                m == j && continue
                term = 1.0
                for k in 1:4
                    (k == j || k == m) || (term *= -nodes[k])
                end
                numerator += term
            end
            weights[j] = numerator / (denominator * scale)
        end

        derivative[i] = sum(weights .* y[indices])
    end
    return derivative
end

function plot_thermodynamics(
    results_path::AbstractString = joinpath(@__DIR__, "results.dat"),
    plot_path::AbstractString = joinpath(@__DIR__, "thermal_properties.png"),
)
    isfile(results_path) || throw(ArgumentError("results file does not exist: $results_path"))
    data = readdlm(results_path, Float64)
    ndims(data) == 2 && size(data, 1) >= 1 ||
        throw(ArgumentError("thermal property plots require at least one row of data: $results_path"))
    size(data, 2) == 3 ||
        throw(ArgumentError("thermal property plots require results columns β F U"))
    all(isfinite, data) ||
        throw(ArgumentError("results file contains non-finite values"))
    all(>(0), data[:, 1]) ||
        throw(ArgumentError("inverse temperatures β must be positive"))

    temperature = 1 ./ data[:, 1]
    order = sortperm(temperature)
    temperature = temperature[order]
    all(diff(temperature) .> 0) ||
        throw(ArgumentError("temperatures must be distinct to calculate thermal properties"))

    free_energy = data[order, 2]
    energy = data[order, 3]
    specific_heat = length(temperature) >= 4 ?
                    four_point_derivative(temperature, energy) : nothing

    free_energy_plot = plot(
        temperature,
        free_energy;
        xlabel = "Temperature T",
        ylabel = "Free energy F",
        label = "F",
        marker = :circle,
        title = "Thermal properties",
        xscale = :log10,
    )
    energy_plot = plot(
        temperature,
        energy;
        xlabel = "Temperature T",
        ylabel = "Internal energy U",
        label = "U",
        marker = :circle,
        xscale = :log10,
    )
    if isnothing(specific_heat)
        figure = plot(
            free_energy_plot,
            energy_plot;
            layout = (2, 1),
            size = (800, 700),
        )
    else
        specific_heat_plot = plot(
            temperature,
            specific_heat;
            xlabel = "Temperature T",
            ylabel = "Specific heat C",
            label = "C = dU/dT",
            marker = :circle,
            xscale = :log10,
        )
        figure = plot(
            free_energy_plot,
            energy_plot,
            specific_heat_plot;
            layout = (3, 1),
            size = (800, 1000),
        )
    end

    savefig(figure, plot_path)
    return (plot = figure, temperature, specific_heat)
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) in (0, 2) ||
        throw(ArgumentError("usage: julia plot_results.jl [results-file plot-file]"))
    plot_thermodynamics(ARGS...)
end
