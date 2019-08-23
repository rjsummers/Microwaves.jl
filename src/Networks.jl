module Networks

using LinearAlgebra

abstract type Network end

mutable struct SParameters <: Network
    """
    S-parameters of a network.
    """
    s::Array{ComplexF64}
    f::Vector{Float64}
    z0::Vector{ComplexF64}
    function SParameters(s, f, z0::Real)
        new(s, f, ones(size(s)[2]).*z0)
    end
end

mutable struct ABCDParameters <: Network
    """
    ABCD-parameters of a network.
    """
    abcd::Array{ComplexF64}
    f::Vector{Float64}
end

mutable struct ZParameters <: Network
    """
    Z-parameters of a network.
    """
    z::Array{ComplexF64}
    f::Vector{Float64}
end

mutable struct YParameters <: Network
    """
    Y-parameters of a network.
    """
    y::Array{ComplexF64}
    f::Vector{Float64}
end

mutable struct HParameters <: Network
    """
    H-parameters of a network.
    """
    h::Array{ComplexF64}
    f::Vector{Float64}
end

mutable struct GParameters <: Network
    """
    G-parameters of a network.
    """
    g::Array{ComplexF64}
    f::Vector{Float64}
end

mutable struct TParameters <: Network
    """
    T-parameters of a network.
    """
    t::Array{ComplexF64}
    f::Vector{Float64}
end

function SParameters(s11, s12, s21, s22, f; z0=50)
    s = cat(cat(s11, s21, dims=2), cat(s12, s22, dims=2), dims=3)
    SParameters(s, f, z0)
end

function ABCDParameters(a, b, c, d, f)
    abcd = cat(cat(a, c, dims=2), cat(b, d, dims=2), dims=3)
    ABCDParameters(abcd, f)
end

function ZParameters(z11, z12, z21, z22, f)
    z = cat(cat(z11, z21, dims=2), cat(z12, z22, dims=2), dims=3)
    ZParameters(z, f)
end

function YParameters(y11, y12, y21, y22, f)
    y = cat(cat(y11, y21, dims=2), cat(y12, y22, dims=2), dims=3)
    YParameters(y, f)
end

function HParameters(h11, h12, h21, h22, f)
    h = cat(cat(h11, h21, dims=2), cat(h12, h22, dims=2), dims=3)
    HParameters(h, f)
end

function GParameters(g11, g12, g21, g22, f)
    g = cat(cat(g11, g21, dims=2), cat(g12, g22, dims=2), dims=3)
    HParameters(g, f)
end

function TParameters(t11, t12, t21, t22, f)
    t = cat(cat(t11, t21, dims=2), cat(t12, t22, dims=2), dims=3)
    TParameters(t, f)
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

function SParameters(z::ZParameters; z0=50)
    s = zeros(z.z)
    z0vec = z0.*Vector(ones(size(z.z)[2]))
    syarr = inv(I * sqrt.(z0vec))
    for i=1:size(z.z)[1]
        s[i,:,:] = (syarr * z.z[i,:,:] * syarr - I)*inv(syarr * z.z[i,:,:] * syarr + I)
    end
    SParameters(s, copy(z.f), z0)
end

function SParameters(y::YParameters; z0=50)
    s = zeros(y.y)
    z0vec = z0.*Vector(ones(size(y.y)[2]))
    szarr = I * sqrt.(z0vec)
    for i=1:size(y.y)[1]
        s[i,:,:] = (I - szarr * y.y[i,:,:] * szarr)*inv(I + szarr * y.y[i,:,:] * szarr)
    end
    SParameters(s, copy(y.f), z0)
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

function ABCDParameters(yparams::YParameters)
    y11 = yparams[:,1,1]
    y12 = yparams[:,1,2]
    y21 = yparams[:,2,1]
    y22 = yparams[:,2,2]
    abcd = zeros(yparams.s)
    abcd[:,1,1] = -y22 ./ y11
    abcd[:,1,2] = -1.0 ./ y21
    abcd[:,2,1] = -(y11.*y22 .- y12.*y21) ./ y21
    abcd[:,2,2] = -y11 ./ y21
    ABCDParameters(abcd, copy(yparams.f))
end

function ZParameters(s::SParameters)
    z = zeros(s.s)
    szarr = I * sqrt.(s.z0)
    for i=1:size(z.z)[1]
        z[i,:,:] = szarr * (I + s.s[i,:,:]) * inv(I - s.s[i,:,:]) * szarr
    end
    ZParameters(z, copy(z.f))
end

function ZParameters(abcd::ABCDParameters)
    a = abcd.abcd[:,1,1]
    b = abcd.abcd[:,1,2]
    c = abcd.abcd[:,2,1]
    d = abcd.abcd[:,2,2]
    z = zeros(abcd.abcd)
    z[:,1,1] = a ./ c
    z[:,1,2] = (a.*d .- b.*c) ./ c
    z[:,2,1] = 1.0 ./ c
    z[:,2,2] = d ./ c
    ZParameters(z, copy(abcd.f))
end

function ZParameters(y::YParameters)
    z = zeros(y.y)
    for i=1:size(z)[2]
        z[i,:,:] = inv(y.y[i,:,:])
    end
    ZParameters(z, copy(y.f))
end

function YParameters(s::SParameters)
    y = zeros(s.s)
    syarr = inv(I * sqrt.(s.z0))
    for i=1:size(y.y)[1]
        y[i,:,:] = syarr * (I - s.s[i,:,:]) * inv(I + s.s[i,:,:]) * syarr
    end
    YParameters(y, copy(y.f))
end

function YParameters(abcd::ABCDParameters)
    a = abcd.abcd[:,1,1]
    b = abcd.abcd[:,1,2]
    c = abcd.abcd[:,2,1]
    d = abcd.abcd[:,2,2]
    y = zeros(abcd.abcd)
    y[:,1,1] = d ./ b
    y[:,1,2] = (b.*c .- a.*d) ./ b
    y[:,2,1] = -1.0 ./ b
    y[:,2,2] =  a ./ b
end

function YParameters(z::ZParameters)
    y = zeros(z.z)
    for i=1:size(z)[2]
        y[i,:,:] = inv(z.z[i,:,:])
    end
    YParameters(y, copy(z.f))
end

end # module