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
end

end # module
