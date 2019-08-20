module RFTools

using SpecialFunctions

function bartlett(N::Int64)
    window = zeros(N)
    for n=0:(N-1)
        window[n+1] = 1 - 2*abs(n - (N-1)/2)/(N-1)
    end
    if N == 1
        window = 1.0
    end
    window
end

function blackman(N::Int64)
    """
    Returns an N-point Blackman window.

    Need to implement the choice of periodic/symmetric Blackman windows.
    """
    if iseven(N)
        M = N / 2
        n = vcat(0:(M-1), (M-1):-1:0)
    else
        M = (N + 1) / 2
        n = vcat(0:(M-1), (M-2):-1:0)
    end
    window = zeros(Float64, N)
    for i=1:N
        window[i] = 0.42 - 0.5.*cos(2*π*n[i]/(N-1)) + 0.08*cos(4*π*n[i]/(N-1))
    end
    window
end

function hamming(N::Int64)
    """
    Returns an N-point hamming window.
    """
    window = zeros(N)
    for n=0:(N-1)
        window[n+1] = 0.54 - 0.46*cos(2*π*n/(N-1))
    end
    window
end

function hanning(N::Int64)
    """
    Returns an N-point hanning window.
    """
    n = 0:(N-1)
    0.5 .* (1 .- cos.(2.0 .* π .* n ./ (N - 1)))
end

function kaiser(N::Int64, β::Float64)
    """
    Returns an N-point Kaiser window with shape factor beta.
    """
    window = zeros(N+1)
    for n=0:N
        window[n+1] = SpecialFunctions.besseli(0, β*sqrt(1 - (2*n/N - 1)^2)) / SpecialFunctions.besseli(0, β)
    end
    window
end

function tukey(N::Int64, α::Float64)
    """
    Returns an N-point Tukey window.
    """
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
