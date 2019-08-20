module Networks

abstract type Network end

mutable struct SParameters <: Network
    s::Array{ComplexF64}
    f::Vector{Float64}
    z0::Array{ComplexF64}
    function SParameters(s, f, z0::Real)
        new(s, f, ones(size(s))*z0)
    end
end

mutable struct ABCDParameters <: Network
    abcd::Array{ComplexF64}
    f::Vector{Float64}
end

mutable struct ZParameters <: Network
    z::Array{ComplexF64}
    f::Vector{Float64}
end

mutable struct YParameters <: Network
    y::Array{ComplexF64}
    f::Vector{Float64}
end

mutable struct HParameters <: Network
    h::Array{ComplexF64}
    f::Vector{Float64}
end

mutable struct GParameters <: Network
    g::Array{ComplexF64}
    f::Vector{Float64}
end

mutable struct TParameters <: Network
    t::Array{ComplexF64}
    f::Vector{Float64}
end

function ABCDParameters(sparams::SParameters)
    s11 = sparams[:,1,1]
    s12 = sparams[:,1,2]
    s21 = sparams[:,2,1]
    s22 = sparams[:,2,2]
    z0 = sparams.z0[1,1,1] # Assume z0 is constant for all ports - this could be extended in the future.
    abcd = zeros(sparams.s)
    abcd[:,1,1] = (1 .+ s11).*(1 .- s22) .+ s12.*s21 ./ (2.0 .* s21)
    abcd[:,1,2] = z0 .* (1 .+ s11).*(1 .+ s22) .- s12.*s21 ./ (2.0 .* s21)
    abcd[:,2,1] = (1/z0) .* (1 .- s11).*(1 .- s22) .- s12.*s21 ./ (2.0 .* s21)
    abcd[:,2,2] = (1 .- s11).*(1 .+ s22) .+ s12.*s21 ./ (2.0 .* s21)
    ABCDParameters(abcd, copy(sparams.f))
end

function ABCDParameters(zparams::ZParameters)
    z11 = zparams[:,1,1]
    z12 = zparams[:,1,2]
    z21 = zparams[:,2,1]
    z22 = zparams[:,2,2]
    abcd = zeros(zparams.s)
    abcd[:,1,1] = z11 ./ z21
    abcd[:,1,2] = (z11.*z22 .- z12.*z21) ./ z21
    abcd[:,2,1] = 1.0 ./ z21
    abcd[:,2,2] = z22 ./ z21
    ABCDParameters(abcd, copy(zparams.f))
end

function SParameters(abcd::ABCDParameters; z0=50)
    a = abcd.abcd[:,1,1]
    b = abcd.abcd[:,1,2]
    c = abcd.abcd[:,2,1]
    d = abcd.abcd[:,2,2]
    s = zeros(abcd.abcd)
    s[:,1,1] = (a .+ b ./ z0 .- c .* z0 .- d) ./ (a .+ b ./ z0 .+ c .* z0 .+ d)
    s[:,1,2] = 2.0 .*(a .* d .- b .* c) ./ (a .+ b ./ z0 .+ c .* z0 .+ d)
    s[:,2,1] = 2.0 ./ (a .+ b ./ z0 .+ c .* z0 .+ d)
    s[:,2,2] = (-a .+ b ./ z0 .- c .* z0 .+ d) ./ (a .+ b ./ z0 .+ c .* z0 .+ d)
    z0arr = ones(size(abcd)) .* z0
    SParameters(s, copy(abcd.f), z0)
end

end # module