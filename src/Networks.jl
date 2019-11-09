module Networks

using LinearAlgebra

abstract type Network end

"""
    SParameters(s, f, z0)

S-parameters of a network.
"""
mutable struct SParameters <: Network
    s::Array{ComplexF64}
    f::Vector{Float64}
    z0::Vector{ComplexF64}
    nports::Int64
    function SParameters(s::Array{<:Complex}, f::Array{<:Real}, z0::Array{<:Real})
        nports = size(s)[2]
        new(s, f, z0, nports)
    end
    function SParameters(s::Array{<:Complex}, f::Array{<:Real}, z0::Real)
        new(s, f, ones(size(s)[2]) .* z0)
    end
end

"""
    ABCDParameters(abcd, f)

ABCD-parameters of a network.
"""
mutable struct ABCDParameters <: Network
    abcd::Array{ComplexF64}
    f::Vector{Float64}
    nports::Int64
    function ABCDParameters(abcd::Array{Complex}, f::Array{Real})
        nports = size(abcd)[2]
        new(abcd, f, nports)
    end
end

"""
    ZParameters(z, f)

Z-parameters of a network.
"""
mutable struct ZParameters <: Network
    z::Array{ComplexF64}
    f::Vector{Float64}
    nports::Int64
    function ZParameters(z::Array{Complex}, f::Array{Real})
        nports = size(z)[2]
        new(z, f, nports)
    end
end

"""
    YParameters(y, f)

Y-parameters of a network.
"""
mutable struct YParameters <: Network
    y::Array{ComplexF64}
    f::Vector{Float64}
    nports::Int64
    function YParameters(y::Array{Complex}, f::Array{Real})
        nports = size(y)[2]
        new(y, f, nports)
    end
end

"""
    HParameters(h, f)

H-parameters of a network.
"""
mutable struct HParameters <: Network
    h::Array{ComplexF64}
    f::Vector{Float64}
    nports::Int64
    function HParameters(h::Array{Complex}, f::Array{Real})
        nports = size(h)[2]
        new(h, f, nports)
    end
end

"""
    GParameters(g, f)

G-parameters of a network.
"""
mutable struct GParameters <: Network
    g::Array{ComplexF64}
    f::Vector{Float64}
    nports::Int64
    function GParameters(g::Array{Complex}, f::Array{Real})
        nports = size(g)[2]
        new(g, f, nports)
    end
end

"""
    TParameters(t, f)

T-parameters of a network.
"""
mutable struct TParameters <: Network
    t::Array{ComplexF64}
    f::Vector{Float64}
    z0::Vector{ComplexF64}
    nports::Int64
    function TParameters(t::Array{Complex}, f::Array{Real}, z0::Array{Real})
        nports = size(t)[2]
        new(t, f, nports)
    end
    function TParameters(t::Array{Complex}, f::Array{Real}, z0::Real)
        new(t, f, ones(size(t)[2]) .* z0)
    end
end

"""
    SParameters(s::SParameters)

Return a copy of the network.

Intended for applications where specific network parameters are desired and
any network type could be provided.
"""
function SParameters(s::SParameters)
    copy(s)
end

"""
    ABCDParameters(abcd::ABCDParameters)

Return a copy of the network.

Intended for applications where specific network parameters are desired and
any network type could be provided.
"""
function ABCDParameters(abcd::ABCDParameters)
    copy(abcd)
end

"""
    ZParameters(z::ZParameters)

Return a copy of the network.

Intended for applications where specific network parameters are desired and
any network type could be provided.
"""
function ZParameters(z::ZParameters)
    copy(z)
end

"""
    YParameters(y::YParameters)

Return a copy of the network.

Intended for applications where specific network parameters are desired and
any network type could be provided.
"""
function YParameters(y::YParameters)
    copy(y)
end

"""
    HParameters(h::HParameters)

Return a copy of the network.

Intended for applications where specific network parameters are desired and
any network type could be provided.
"""
function HParameters(h::HParameters)
    copy(h)
end

"""
    GParameters(g::GParameters)

Return a copy of the network.

Intended for applications where specific network parameters are desired and
any network type could be provided.
"""
function GParameters(g::GParameters)
    copy(g)
end

"""
    TParameters(t::TParameters)

Return a copy of the network.

Intended for applications where specific network parameters are desired and
any network type could be provided.
"""
function TParameters(t::TParameters)
    copy(t)
end

"""
    SParameters(s11, s12, s21, s22, f[, z0=50])

Creates a two-port S-parameter network given four vectors
for each parameter vs. frequency.
"""
function SParameters(s11, s12, s21, s22, f; z0 = 50)
    s = cat(cat(s11, s21, dims = 2), cat(s12, s22, dims = 2), dims = 3)
    SParameters(s, f, z0)
end

"""
    ABCDParameters(a, b, c, d, f)

Creates a two-port ABCD-parameter network given four vectors
for each parameter vs. frequency.
"""
function ABCDParameters(a, b, c, d, f)
    abcd = cat(cat(a, c, dims = 2), cat(b, d, dims = 2), dims = 3)
    ABCDParameters(abcd, f)
end

"""
    ZParameters(z11, z12, z21, z22, f)

Creates a two-port Z-parameter network given four vectors
for each parameter vs. frequency.
"""
function ZParameters(z11, z12, z21, z22, f)
    z = cat(cat(z11, z21, dims = 2), cat(z12, z22, dims = 2), dims = 3)
    ZParameters(z, f)
end

"""
    YParameters(y11, y12, y21, y22, f)

Creates a two-port Y-parameter network given four vectors
for each parameter vs. frequency.
"""
function YParameters(y11, y12, y21, y22, f)
    y = cat(cat(y11, y21, dims = 2), cat(y12, y22, dims = 2), dims = 3)
    YParameters(y, f)
end

"""
    HParameters(h11, h12, h21, h22, f)

Creates a two-port H-parameter network given four vectors
for each parameter vs. frequency.
"""
function HParameters(h11, h12, h21, h22, f)
    h = cat(cat(h11, h21, dims = 2), cat(h12, h22, dims = 2), dims = 3)
    HParameters(h, f)
end

"""
    GParameters(g11, g12, g21, g22, f)

Creates a two-port G-parameter network given four vectors
for each parameter vs. frequency.
"""
function GParameters(g11, g12, g21, g22, f)
    g = cat(cat(g11, g21, dims = 2), cat(g12, g22, dims = 2), dims = 3)
    HParameters(g, f)
end

"""
    TParameters(t11, t12, t21, t22, f[, z0=50])

Creates a two-port T-parameter network given four vectors
for each parameter vs. frequency.
"""
function TParameters(t11, t12, t21, t22, f; z0 = 50)
    t = cat(cat(t11, t21, dims = 2), cat(t12, t22, dims = 2), dims = 3)
    TParameters(t, f, z0)
end

"""
    SParameters(abcd::ABCDParameters[, z0=50])

Converts an ABCD-parameter object to an S-parameter object.
Only works for two-port networks, and assumes z0 is the same
for all ports.
"""
function SParameters(abcd::ABCDParameters; z0 = 50)
    a = abcd.abcd[:,1,1]
    b = abcd.abcd[:,1,2]
    c = abcd.abcd[:,2,1]
    d = abcd.abcd[:,2,2]
    s = zeros(abcd.abcd)
    s[:,1,1] = (a .+ b ./ z0 .- c .* z0 .- d) ./ (a .+ b ./ z0 .+ c .* z0 .+ d)
    s[:,1,2] = 2.0 .* (a .* d .- b .* c) ./ (a .+ b ./ z0 .+ c .* z0 .+ d)
    s[:,2,1] = 2.0 ./ (a .+ b ./ z0 .+ c .* z0 .+ d)
    s[:,2,2] = (-a .+ b ./ z0 .- c .* z0 .+ d) ./ (a .+ b ./ z0 .+ c .* z0 .+ d)
    z0arr = ones(size(abcd)) .* z0
    SParameters(s, copy(abcd.f), z0)
end

"""
    SParameters(z::ZParameters[, z0=50])

Converts a Z-Parameter object to an S-Parameter object.
"""
function SParameters(z::ZParameters; z0 = 50)
    s = zeros(z.z)
    z0vec = z0 .* Vector(ones(size(z.z)[2]))
    syarr = inv(I * sqrt.(z0vec))
    for i = 1:size(z.z)[1]
        s[i,:,:] = (syarr * z.z[i,:,:] * syarr - I) * inv(syarr * z.z[i,:,:] * syarr + I)
    end
    SParameters(s, copy(z.f), z0)
end

"""
    SParameters(y::YParameters[, z0=50])

Converts a Y-Parameter object to an S-Parameter object.
"""
function SParameters(y::YParameters; z0 = 50)
    s = zeros(y.y)
    z0vec = z0 .* Vector(ones(size(y.y)[2]))
    szarr = I * sqrt.(z0vec)
    for i = 1:size(y.y)[1]
        s[i,:,:] = (I - szarr * y.y[i,:,:] * szarr) * inv(I + szarr * y.y[i,:,:] * szarr)
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
    abcd[:,1,1] = (1 .+ s11) .* (1 .- s22) .+ s12 .* s21 ./ (2.0 .* s21)
    abcd[:,1,2] = z0 .* (1 .+ s11) .* (1 .+ s22) .- s12 .* s21 ./ (2.0 .* s21)
    abcd[:,2,1] = (1 / z0) .* (1 .- s11) .* (1 .- s22) .- s12 .* s21 ./ (2.0 .* s21)
    abcd[:,2,2] = (1 .- s11) .* (1 .+ s22) .+ s12 .* s21 ./ (2.0 .* s21)
    ABCDParameters(abcd, copy(sparams.f))
end

function ABCDParameters(zparams::ZParameters)
    z11 = zparams[:,1,1]
    z12 = zparams[:,1,2]
    z21 = zparams[:,2,1]
    z22 = zparams[:,2,2]
    abcd = zeros(zparams.s)
    abcd[:,1,1] = z11 ./ z21
    abcd[:,1,2] = (z11 .* z22 .- z12 .* z21) ./ z21
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
    abcd[:,2,1] = -(y11 .* y22 .- y12 .* y21) ./ y21
    abcd[:,2,2] = -y11 ./ y21
    ABCDParameters(abcd, copy(yparams.f))
end

function ZParameters(s::SParameters)
    z = zeros(s.s)
    szarr = I * sqrt.(s.z0)
    for i = 1:size(z.z)[1]
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
    z[:,1,2] = (a .* d .- b .* c) ./ c
    z[:,2,1] = 1.0 ./ c
    z[:,2,2] = d ./ c
    ZParameters(z, copy(abcd.f))
end

function ZParameters(y::YParameters)
    z = zeros(y.y)
    for i = 1:size(z)[2]
        z[i,:,:] = inv(y.y[i,:,:])
    end
    ZParameters(z, copy(y.f))
end

function YParameters(s::SParameters)
    y = zeros(s.s)
    syarr = inv(I * sqrt.(s.z0))
    for i = 1:size(y.y)[1]
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
    y[:,1,2] = (b .* c .- a .* d) ./ b
    y[:,2,1] = -1.0 ./ b
    y[:,2,2] =  a ./ b
end

function YParameters(z::ZParameters)
    y = zeros(z.z)
    for i = 1:size(z)[2]
        y[i,:,:] = inv(z.z[i,:,:])
    end
    YParameters(y, copy(z.f))
end

function combine_touchstone_nums(x::Real, y::Real, format::String)
    if format == "DB"
        return (10^(x / 20)) * exp(1im * y * π / 180)
    elseif format == "MA"
        return x * exp(1im * y * π / 180)
    elseif format == "RI"
        return complex(x, y)
    end
end

function parse_touchstone_opt_line(row::String)
    funit = 1e9
    parameter = "S"
    format = "MA"
    z0 = 50.0
    next_option_is_z0 = false

    for option = split(row)[2:end]
        if option == "Hz"
            funit = 1e0
        elseif option == "kHz"
            funit = 1e3
        elseif option == "MHz"
            funit = 1e6
        elseif option == "GHz"
            funit = 1e9
        elseif occursin(option, "SZYHG")
            parameter = option
        elseif option == "DB"
            format = string(option)
        elseif option == "MA"
            format = string(option)
        elseif option == "RI"
            format = string(option)
        elseif option == "R"
            next_option_is_z0 = true
        elseif next_option_is_z0
            z0 = parse(Float64, option)
            next_option_is_z0 = false
        end
    end
    (funit, parameter, format, z0)
end

"""
    read_touchstone(filename::String)

Reads a touchstone file and returns a corresponding network object.

Currently swaps axis order used in rest of Networks.jl. I plan to update the
rest of Networks.jl soon to match the convention used in MatLab's RF toolbox.

Works on version 1 files with number of ports <= 2. Larger numbers of ports
might work, but might also have jumbled up parameters.
"""
function read_touchstone(filename::String)
    m = match(r"(\d+)(?!.*\d)", filename)

        freq = zeros(0)
        data = zeros(0)
        opts = ()
        open(filename) do file
            for ln in eachline(file)
        # Apparently skips comments and whitespace
        # Taken from https://github.com/JuliaIO/ConfParser.jl/blob/master/src/ConfParser.jl
                occursin(r"^\s*(\n|\!|;)", ln) && continue
                occursin(r"\w", ln) || continue

                if ln[1] == '#'
                    opts = parse_touchstone_opt_line(ln)
                    continue
                end

                rowdata = parse.(Float64, split(ln))
                if  isodd(length(rowdata))
                    append!(freq, rowdata[1])
                    append!(data, rowdata[2:end])
                else
                    append!(data, rowdata)
                end
            end
        end
        z0 = opts[4]
        parameter = opts[2]
        nports = div(length(data), 4*length(freq))
        params = combine_touchstone_nums.(data[1:2:end], data[2:2:end], opts[3])
        params = reshape(params, nports, nports, :)



    if parameter == "S"
        net = SParameters(params, freq, z0)
    elseif parameter == "Z"
        net = ZParameters(params, freq)
    elseif parameter == "Y"
        net = YParameters(params, freq)
    elseif parameter == "H"
        net = HParameters(params, freq)
    elseif parameter == "G"
        net = GParameters(params, freq)
    end
    net
end

end # module
