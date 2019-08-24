__precompile__(true)

module RFTools

include("Networks.jl")
include("Windows.jl")

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

function voltage_gain(x::Number, units::Symbol)
    if units == :db
        return 10^(x/20)
    elseif units == :np
        return exp(x)
    end
end

end # module
