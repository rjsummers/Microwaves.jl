__precompile__(true)

module RFTools

include("Networks.jl")
include("Windows.jl")

function voltage_gain(x::Number, units::Symbol)
    if units == :power
        return sqrt(x)
    elseif units == :db
        return 10^(x/20)
    elseif units == :np
        return exp(x)
    end
end

function power_gain(x::Number, units::Symbol)
    if units == :voltage
        return x^2
    elseif units == :db
        return 10^(x/10)
    elseif units == :np
        return 10^(20 * log10(Base.MathConstants.e) * x / 10)
    end
end

function db(x::Number)
    20*log10(abs(x))
end

function db(x::Number, units::Symbol)
    if units == :voltage
        return 20*log10(abs(x))
    elseif units == :power
        return 10*log10(x)
    elseif units == :np
        return 20 * log10(Base.MathConstants.e) * x
    end
end 

end # module
