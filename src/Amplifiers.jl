module Amplifiers

using ..RF: Networks

function Γ_in(net::Network, Γ_L::Vector{Complex})
    sparams = SParameters(net)
    sparams[:,1,1] .+ (sparams[:,1,2] .* sparams[:,2,1] .* Γ_L) ./ (1 .- sparams[:,2,2] .* Γ_L)
end

function Γ_out(net::Network, Γ_S::Vector{Complex})
    sparams = SParameters(net)
    sparams[:,2,2] .+ (sparams[:,1,2] .* sparams[:,2,1] .* Γ_S) ./ (1 .- sparams[:,1,1] .* Γ_S)
end

end # module
