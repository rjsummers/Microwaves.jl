module Amplifiers

using ..RF.Networks: Network, SParameters

function Γ_in(net::Network, Γ_L::Vector{Complex})
    sparams = SParameters(net)
    sparams.s[:,1,1] .+ (sparams.s[:,1,2] .* sparams.s[:,2,1] .* Γ_L) ./ (1 .- sparams.s[:,2,2] .* Γ_L)
end

function Γ_out(net::Network, Γ_S::Vector{Complex})
    sparams = SParameters(net)
    sparams.s[:,2,2] .+ (sparams.s[:,1,2] .* sparams.s[:,2,1] .* Γ_S) ./ (1 .- sparams.s[:,1,1] .* Γ_S)
end

function power_gain(net::Network, Γ_L::Vector{Complex})
    s = SParameters(net)
    num = abs.(s.s[:,2,1]).^2 .* (1 .- abs.(Γ_L).^2)
    den = (1 .- abs.(Γ_in(s, Γ_L)).^2) .* abs.(1 .- s.s[:,2,2] .* Γ_L).^2
    num ./ den
end

function available_power_gain(net::Network, Γ_S::Vector{Complex})
    s = SParameters(net)
    num = abs.(s.s[:,2,1]).^2 .* (1 .- abs.(Γ_S).^2)
    den = abs.(1 .- s.s[:,1,1] .* Γ_s).^2 .* (1 .- Γ_out(s, Γ_S).^2)
    num ./ den
end

function transducer_power_gain(net::Network, Γ_S::Vector{Complex}, Γ_L::Vector{Complex})
    s = SParameters(net)
    num = abs.(s.s[:,2,1]).^2 .* (1 .- abs.(Γ_S).^2) .* (1 .- abs.(Γ_L).^2)
    den = abs.(1 .- Γ_S .* Γ_in(s, Γ_L)).^2 .* abs(1 .- s.s[:,2,2] .* Γ_L).^2
    num ./ den
end

function μ_test(net::Network; retmu=false)
    s = SParameters(net)
    Δ = s.s[:,1,1] .* s.s[:,2,2] .- s.s[:,1,2] .* s.s[:,2,1]
    μ = (1 .- abs.(s.s[:,1,1]).^2) ./ (abs.(s.s[:,2,2] .- Δ .* conj.(s.s[:,1,1])) .+
                                       abs.(s.s[:,1,2] .* s.s[:,2,1]))
    result = true
    for muval=μ
        if muval <= 1
            result = false
        end
    end
    if retmu
        return result, μ
    end
    result
end

end # module
