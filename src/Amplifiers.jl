module Amplifiers

using ..RF: Networks

function Γ_in(net::Network, Γ_L::Vector{Complex})
    sparams = SParameters(net)
    sparams[:,1,1] .+ (sparams[:,1,2] .* sparams[:,2,1] .* Γ_L) ./ (1 .- sparams[:,2,2] .* Γ_L)
end

end # module