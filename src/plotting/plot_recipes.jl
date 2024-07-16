using RecipesBase

DEFAULT_ALPHA = 0.2
DEFAULT_LABEL = ""
DEFAULT_GRID = false
DEFAULT_LEGEND = false
DEFAULT_DISTRIBUTION = :pdf
DEFAULT_FILL = :gray
DEFAULT_COLOUR_PDF = :blue
DEFAULT_COLOUR_UPPER = :red
DEFAULT_COLOUR_LOWER = :black

DEFAULT_PLOT_RANGE_EXTEND = 0.2
DEFAULT_PLOT_RANGE_INTERVAL = 0.4
DEFAULT_PLOT_GRID_NUMBER = 500
DEFAULT_FONT_SIZE = 18
DEFAULT_TICK_SIZE = 12


@recipe function _plot(x::RandomVariable)

    grid --> DEFAULT_GRID
    legend --> DEFAULT_LEGEND
    ylabel --> "pdf"
    xlabel --> x.name

    lo_grid = quantile(x, 0.001)
    hi_grid = quantile(x, 0.999)

    lo_grid = x.lb
    hi_grid = x.ub

    width = (hi_grid + lo_grid)/2

    lo_grid = lo_grid - abs(width * DEFAULT_PLOT_RANGE_EXTEND)
    hi_grid = hi_grid + abs(width * DEFAULT_PLOT_RANGE_EXTEND)

    x_grid = range(lo_grid, hi_grid, DEFAULT_PLOT_GRID_NUMBER)
    pdf_evals = pdf.(Ref(x), x_grid)
    
    @series begin
        fillrange := pdf_evals
        color := DEFAULT_FILL
        fillalpha := DEFAULT_ALPHA
        x_grid, zeros(length(x_grid))
    end

    @series begin
        color --> DEFAULT_COLOUR_PDF
        alpha := 1
        label := ""
        x_grid, pdf_evals
    end
end

@recipe function _plot(x::Interval)

    grid --> DEFAULT_GRID
    legend --> DEFAULT_LEGEND
    ylabel --> "cdf"
    xlabel --> x.name

    lo_grid = x.lb
    hi_grid = x.ub

    width = (hi_grid + lo_grid)/2

    lo_grid = lo_grid - abs(width * DEFAULT_PLOT_RANGE_INTERVAL)
    hi_grid = hi_grid + abs(width * DEFAULT_PLOT_RANGE_INTERVAL)

    x_grid = range(lo_grid, hi_grid, DEFAULT_PLOT_GRID_NUMBER)

    cdf_lo = x_grid .>= x.ub
    cdf_hi = x_grid .>= x.lb
    
    @series begin
        fillrange := cdf_hi
        color := DEFAULT_FILL
        fillalpha := DEFAULT_ALPHA
        x_grid, cdf_lo
    end

    @series begin
        color --> DEFAULT_COLOUR_LOWER
        alpha := 1
        label := ""
        x_grid, cdf_lo
    end

    @series begin
        color --> DEFAULT_COLOUR_UPPER
        alpha := 1
        label := ""
        x_grid, cdf_hi
    end

end

@recipe function _plot(x::ProbabilityBox)

    grid --> DEFAULT_GRID
    legend --> DEFAULT_LEGEND
    ylabel --> "cdf"
    xlabel --> x.name

    lo_grid = quantile(x, 0.001)[1]
    hi_grid = quantile(x, 0.999)[2]

    width = (hi_grid + lo_grid)/2

    lo_grid = lo_grid - abs(width * DEFAULT_PLOT_RANGE_EXTEND)
    hi_grid = hi_grid + abs(width * DEFAULT_PLOT_RANGE_EXTEND)

    x_grid = range(lo_grid, hi_grid, DEFAULT_PLOT_GRID_NUMBER)
    cdf_evals = cdf.(Ref(x), x_grid)

    @series begin
        fillrange := hi.(cdf_evals)
        color := DEFAULT_FILL
        fillalpha := DEFAULT_ALPHA
        x_grid, lo.(cdf_evals)
    end

    @series begin
        color --> DEFAULT_COLOUR_LOWER
        alpha := 1
        label := ""
        x_grid, lo.(cdf_evals)
    end

    @series begin
        color --> DEFAULT_COLOUR_UPPER
        alpha := 1
        label := ""
        x_grid, hi.(cdf_evals)
    end
end


###
#   This code is a modified version of the plot recipe from IntervalArithmetic.jl
#       https://github.com/JuliaIntervals/IntervalArithmetic.jl       
###

# Plot a 2D IntervalBox:
@recipe function _plot(x::Interval, y::Interval) #; customcolor = :green, alpha=0.5)

    seriesalpha --> DEFAULT_ALPHA
    seriestype := :shape

    x = [x.lb, x.ub, x.ub, x.lb]
    y = [y.lb, y.lb, y.ub, y.ub]

    x, y
end

# Plot a vector of 2D IntervalBoxes:
@recipe function _plot(xx::Vector{T}, yy::Vector{T}) where T<:Interval

    seriestype := :shape

    xs = Float64[]
    ys = Float64[]

    # build up coords:  # (alternative: use @series)
    for i = 1:length(xx)
        (x, y) = (xx[i], yy[i])

        # use NaNs to separate
        append!(xs, [x.lb, x.ub, x.ub, x.lb, NaN])
        append!(ys, [y.lb, y.lb, y.ub, y.ub, NaN])

    end

    seriesalpha --> DEFAULT_ALPHA
    xs, ys

end
