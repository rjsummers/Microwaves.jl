module RFTools

function tukey(N::Int64, α::Float64)
    window = zeros(Float64, N)
    boundary1 = α*N/2
    boundary2 = N*(1-α/2)
    for n=0:(N-1)
        if n < boundary1
            window[n+1] = 0.5*(1 + cos(π*(2*n/(α*N) - 1)))
        elseif n <= boundary2
            window[n+1] = 1
        else
            window[n+1] = 0.5*(1 + cos(π*(2*n/(α*N) - 2/α + 1)))
        end
    end
    window
end

end # module
