module Networks

using LinearAlgebra

abstract type Network end

"""
S-parameters of a network.
"""
mutable struct SParameters <: Network
    s::Array{ComplexF64}
    f::Vector{Float64}
    z0::Vector{ComplexF64}
    function SParameters(s, f, z0::Real)
        new(s, f, ones(size(s)[2]).*z0)
    end
end

"""
ABCD-parameters of a network.
"""
mutable struct ABCDParameters <: Network
    abcd::Array{ComplexF64}
    f::Vector{Float64}
end

"""
Z-parameters of a network.
"""
mutable struct ZParameters <: Network
    z::Array{ComplexF64}
    f::Vector{Float64}
end

"""
Y-parameters of a network.
"""
mutable struct YParameters <: Network
    y::Array{ComplexF64}
    f::Vector{Float64}
end

"""
H-parameters of a network.
"""
mutable struct HParameters <: Network
    h::Array{ComplexF64}
    f::Vector{Float64}
end

"""
G-parameters of a network.
"""
mutable struct GParameters <: Network
    g::Array{ComplexF64}
    f::Vector{Float64}
end

"""
T-parameters of a network.
"""
mutable struct TParameters <: Network
    t::Array{ComplexF64}
    f::Vector{Float64}
    z0::Vector{ComplexF64}
    function SParameters(s, f, z0::Real)
        new(s, f, ones(size(s)[2]).*z0)
    end
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

function TParameters(t11, t12, t21, t22, f; z0=50)
    t = cat(cat(t11, t21, dims=2), cat(t12, t22, dims=2), dims=3)
    TParameters(t, f, z0)
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

"""
    read_touchstone(filename::String)

Reads a touchstone file and returns a corresponding network object.

Currently, only works on 1 and 2 port version 1 files,
with a broadening of capability coming in the near future.
"""
function read_touchstone(filename::String)
    function combine_nums(x::Real, y::Real, format::String)
        if format == "DB"
            return (10^(x/20)) * exp(1im * y * π / 180)
        elseif format == "MA"
            return x * exp(1im * y * π / 180)
        elseif format == "RI"
            return complex(x, y)
        end
    end

    function parse_options(row::String)
        funit = 1e9
        parameter = "S"
        format = "MA"
        z0 = 50.0
        next_option_is_z0 = false

        for option=split(row, " ")[2:end]
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

    contents = open("mwo_butter.s2p") do file
        [line for line in eachline(file)]
    end

    # Remove comments and empty lines
    nlines = length(contents)
    for i=1:nlines
        contents[i] = strip(split(contents[i], "!")[1])
    end
    filter!(line -> !isempty(line), contents)

    if contents[1][1] == '#'
        version = 1.0
        (funit, parameter, format, z0) = parse_options(contents[1])
        data_after = 1
        freq = Vector{Float64}(undef, 0)
        data_vector = Vector{ComplexF64}(undef, 0)
    elseif occursin("[Version]", contents[1])
        version = 2.0
        (funit, parameter, format, z0) = parse_options(contents[2])

        # Gather data from v2 keywords
        reading_reference = false
        refs_read = 0
        mformat = "Full"
        for i=3:length(contents)
            row = contents[i]
            if occursin("[Number of Ports]", row)
                nports = parse(Int64, split(row, " ")[end])
            elseif occursin("[Two-Port Data Order]", row) && (nports == 2)
                two_port_order = split(row, " ")[end]
            elseif occursin("[Number of Frequencies]", row)
                nfreq = parse(Int64, split(row, " ")[end])
            elseif occursin("[Number of Noise Frequencies]", row)
                nnoisefreq = parse(Int64, split(row, " ")[end])
            elseif occursin("[Reference]", row)
                reading_reference = true
                refs_read = 0
                z0 = zeros(nports)
                refs = split(row, " ")[2:end]
                for i=2:length(refs)
                    z0[i] = parse(Float64, refs[i])
                    refs_read += 1
                end
                if refs_read == nports
                    reading_reference = false
                end
            elseif occursin("[Matrix Format]", row)
                mformat = split(row, " ")[end]
            elseif occursin("[Network Data]", row)
                data_after = i
                freq = Vector{Float64}(undef, nfreq)
                data_vector = Vector{ComplexF64}(undef, nfreq * nports^2)
            elseif reading_reference
                refs = split(row, " ")
                for i=1:length(refs)
                    z0[refs_read + 1] = parse(Float64, refs[i])
                    refs_read == 1
                end
                if refs_read == nports
                    reading_reference = false
                end
            end
        end
    end

    if version == 1.0
        first_data_row = true
        freq_row_length = -1

        for row=contents[(data_after + 1):end]
            data = parse.(Float64, split(row, " "))

            if first_data_row
                freq_row_length = length(data)
                first_data_row = false
            end

            if length(data) == freq_row_length
                push!(freq, data[1] * funit)
                for i=2:2:(length(data) - 1)
                    push!(data_vector, combine_nums(data[i], data[i + 1], format))
                end
            else
                for i=1:2:(length(data) - 1)
                    push!(data_vector, combine_nums(data[i], data[i + 1], format))
                end
            end
        end
    end

    if version == 1.0
        nfreq = length(freq)
        nports = floor(Int64, sqrt(length(data_vector) / nfreq))
        if nports^2 * nfreq != length(data_vector)
            error("Data points do not match frequencies.")
        end
        if nports == 1
            data = reshape(data_vector, nfreq, 1, 1)
        elseif nports == 2
            data = zeros(ComplexF64, nfreq, 2, 2)
            for i=1:nfreq
                data[i,1,1] = data_vector[4*(i-1) + 1]
                data[i,2,1] = data_vector[4*(i-1) + 2]
                data[i,1,2] = data_vector[4*(i-1) + 3]
                data[i,2,2] = data_vector[4*(i-1) + 4]
            end
        end
    end

    if parameter == "S"
        net = SParameters(data, freq, z0)
    elseif parameter == "Z"
        net = ZParameters(data, freq)
    elseif parameter == "Y"
        net = YParameters(data, freq)
    elseif parameter == "H"
        net = HParameters(data, freq)
    elseif parameter == "G"
        net = GParameters(data, freq)
    end
    net
end

end # module
