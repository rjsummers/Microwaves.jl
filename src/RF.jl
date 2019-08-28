__precompile__(true)

module RF

include("Networks.jl")

function voltage_gain(x::Number, units::Symbol)
    if units == :power
        return sqrt(x)
    elseif units == :dB
        return 10^(x/20)
    elseif units == :Np
        return exp(x)
    end
end

function power_gain(x::Number, units::Symbol)
    if units == :voltage
        return x^2
    elseif units == :dB
        return 10^(x/10)
    elseif units == :Np
        return 10^(20 * log10(Base.MathConstants.e) * x / 10)
    end
end

function dB(x::Number)
    20*log10(abs(x))
end

function dB(x::Number, units::Symbol)
    if units == :voltage
        return 20*log10(abs(x))
    elseif units == :power
        return 10*log10(x)
    elseif units == :Np
        return 20 * log10(Base.MathConstants.e) * x
    end
end 

function Np(x::Number, units::Symbol)
    if units == :voltage
        return log(x)
    elseif units == :power
        return log(sqrt(x))
    elseif units == :dB
        return (1.0/20.0) * log(10) * x
    end
end

end # module
