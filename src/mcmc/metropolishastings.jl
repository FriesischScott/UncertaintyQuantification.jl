function metropolishastings(U0::AbstractMatrix, proposal::Sampleable{Univariate})
    dim = size(U0)
    Φ = MvNormal(dim[2], 1.0)

    U = rand(proposal, dim...)

    pdf_i = prod(pdf(Φ, (U0 .+ U)'), dims=2)
    pdf_0 = prod(pdf(Φ, U0'), dims=2)

    α = pdf_i ./ pdf_0
    accept = α .>= rand(size(α)...)

    U0 += U .* accept
end