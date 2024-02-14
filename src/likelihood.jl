# using SpecialFunctions
# function poisson_nloglikelihood(k, λ) # element by element, evaluate -log Lp when k is expected to be =poisson(λ), i.e. k=data, λ=model
# return λ .-k.*log.(λ)+lgamma.(k .+ 1.0)
# end

using Distributions
function poisson_likelihood_image(x, y)
return -sum(loglikelihood.(Poisson.(y), x))
end

function add_poisson_noise!(x::Vector{Float64})
  x[:] += rand.(Poisson.(x))    
end