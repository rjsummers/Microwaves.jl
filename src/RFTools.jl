module RFTools

function tukey(N::Int64, α::Float64)
    window = zeros(Float64, N)
    boundary1 = α/2
    boundary2 = 1 - α/2
    x=range(0, length=N, stop=1)
    for n=1:N
        if x[n] < boundary1
            window[n] = 0.5*(1 + cos(2*π/α*(x[n] - α/2)))
        elseif x[n] < boundary2
            window[n] = 1
        else
            window[n] = 0.5*(1 + cos(2*π/α*(x[n] - 1 + α/2)))
        end
    end
    window
end

end # module
