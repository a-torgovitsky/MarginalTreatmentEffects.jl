# TODO(omkarakatta): consider getting rid of `isnothing(steps)` conditional?

"""
    df_to_coordinates(df::DataFrame,
                      xindex,
                      yindex;
                      steps::Union{Number, Nothing} = nothing,
                      disctol::Number = 7)

This function returns a vector of strings. Each string is of the form
    "(x, y)(x, y) ⋯ (x, y)"
where each (x, y) tuple is a coordinate for the `addplot` command.
Each string of coordinates is associated with its own `addplot` command.
Consequently, all the (x, y) tuples within a string will be connected;
the (x, y) tuples across different strings will not be connected.
This design is useful when there are discontinuities.

The primary use of this function is to plot the following:
    1. weights associated with IV-like estimands and target parameters
    2. MTRs

The MTRs and the weights use different strategies for creating their
coordinates because we know the weights are piecewise-constant. Hence, we can
draw them without any sort of approximation. To ensure there are tick marks
within each piece, new points need to be introduced within each piece (see the
`steps` argument).

On the other hand, we don't know where the MTRs are discontinuous.
So, we discover where they have discontinuities by checking the slopes of the
secant lines between two consecutive points. If the magnitude of the slopes is
larger than `disctol`, then we have found a discontinuity. These
discontinuities determine the segments.

# Parameters

- df,xindex,yindex: df[:, xindex] and df[:, yindex] contains (x, y) pairs
- steps: By default, this option is not specified, in which case, no new points
    will be introduced. If a number is specified, then new points will be
    introduced on the line segment between two points:
        (df[i, xindex], df[i, yindex]) and (df[i + 1, xindex], df[i, yindex]).
    The distance between these new points will be `steps` units.
- disctol: If the secant line between two consecutive points is larger than
    this value, then we consider there to be a discontinuity. As a result, a
    new entry will be made in the output vector.
"""
function df_to_coordinates(df::DataFrame,
                           xindex,
                           yindex;
                           steps::Union{Number, Nothing} = nothing,
                           disctol::Number = 7)
    x = round.(df[:, xindex], digits = 4)
    y = round.(df[:, yindex], digits = 4)
    segments = Vector{String}()
    # If `steps` is not specified, return a vector of coordinates, where
    #   discontinuities separate vectors from one another.
    #   This is useful for drawing the MTRs.
    # If `steps` is specified, return a vector of coordinates (i.e., endpoints)
    #   `steps` controls how far apart the points are.
    #   This is useful for drawing piecewise constant functions for the weights.
    if isnothing(steps)
        # Find indices after which discontinuities occur
        lagdiff = v -> v - vcat(v[2:length(v)], 0)
        slope = lagdiff(y) ./ lagdiff(x)
        disc = findall(abs.(slope) .> disctol)
        coordinates = "" # initialize string of coordinates
        for i in 1:nrow(df)
            coordinates = coordinates *
                "(" * string(x[i]) * "," * string(y[i]) * ")"
            # After arriving at a discontinuity, start new string of coordinates
            if i in disc
                push!(segments, coordinates)
                coordinates = ""
            end
            if i == nrow(df)
                push!(segments, coordinates)
            end
        end
    else
        for i in 1:nrow(df)
            yval = y[i]
            if yval ≈ 0 # don't plot trivial values
                continue
            end
            coordinates = ""
            left = x[i]
            if i == nrow(df)
                right = 1
            else
                right = x[i + 1]
            end
            grid = round.(range(left, right, step = steps), digits = 3)
            unique!(push!(grid, right)) # ensure that upper is in grid
            for point_idx in 1:length(grid)
                coordinates = coordinates * "(" *
                    string(grid[point_idx]) * "," *
                    string(yval) * ")"
            end
            push!(segments, coordinates)
        end
    end
    return segments
end

# Generate the title used in the legend
function legendtitle(tp::TargetParameter)
    return tp.legendtitle
end
function legendtitle(ivlike::IVLike)
    return ivlike.legendtitle[1]
end

# Generate the title used in TikZ path for cross-referencing
function pathtitle(tp::TargetParameter)
    if tp.name == "LATE(u₁, u₂)"
        title = "late"
    elseif tp.name == "ATT"
        title = "att"
    else
        @error "Path title for target parameter is not supported." tp.name
    end
    return title
end

function pathtitle(ivlike::IVLike)
    if ivlike.name == "IV Slope"
        title = "ivs"
    elseif ivlike.name == "OLS Slope"
        title = "olss"
    elseif occursin("IV Slope for 𝟙(Z == z) for z ∈", ivlike.name)
        title = "ivnps" .* string.(ivlike.params[:support])
    elseif occursin("TSLS Slope for 𝟙(Z == z) for z ∈", ivlike.name)
        title = "tslss"
    elseif occursin("Wald", ivlike.name)
        title = "wald"
    elseif ivlike.name == "Saturated"
        title = "saturated" .* string.(collect(1:length(ivlike.s)))
    else
        @error "Path title for IV-like estimand is not supported." tp.name
    end
    return title
end

# Rewrite lower and upper bounds as an interval
function parse_bounds(result)
    lb = result[:lb]
    ub = result[:ub]
    return ": [\$ $(@sprintf("%.3f", lb)), $(@sprintf("%.3f", ub)) \$]"
end
