function copularand(covariance::Matrix{<:Number}, n::Int64, m::Int64)
    a = cholesky(covariance).L
    z = transpose(rand(Normal(), n, m))
    x = a * z
    u = transpose(cdf.(Normal(), x))
end