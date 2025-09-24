function linear_binning(x::AbstractVector, nbins::Integer)
    xmin = minimum(x)
    xmax = maximum(x)

    g = collect(range(xmin, xmax; length=nbins)) # grid points
    w = zeros(eltype(x), nbins) # weights
    Δ = (xmax - xmin) / (nbins - 1)
    for xi in x
        # Find the g index on the left
        idx = Int(clamp(floor((xi - xmin) / Δ) + 1, 1, nbins))
        if idx == nbins
            # If at the last bin, assign all weight there
            w[nbins] += 1.0
        else
            left = g[idx]
            # Linear weights
            w_right = (xi - left) / Δ
            w_left = 1.0 - w_right
            w[idx] += w_left
            w[idx + 1] += w_right
        end
    end
    return g, w
end
