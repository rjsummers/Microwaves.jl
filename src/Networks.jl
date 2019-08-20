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

# function SParameters(abcd::ABCDParameters; z0=50)
#     a = abcd(:,1,1)
#     b = abcd(:,1,2)
#     c = abcd(:,2,1)
#     d = abcd(:,2,2)
#     s = zeros(abcd)
#     s(:,1,1) = (a + b / z0 - c * z0 - d) / (a + b / z0 + c * z0 + d)
#     s(:,1,2) = 2*(a * d - b * c) / (a + b / z0 + c * z0 + d)
#     s(:,2,1) = 2 / (a + b / z0 + c * z0 + d)
#     s(:,2,2) = (-a + b / z0 - c * z0 + d) / (a + b / z0 + c * z0 + d)
#     z0arr = ones(size(abcd)) * z0
#     SParameters(s, copy(abcd.f), z0)
# end

end # module