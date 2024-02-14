# Note: could use PoissonRandom for Poisson noise instead of Distributions

using Distributions
function poisson_likelihood_image(model, data)
return -sum(loglikelihood.(Poisson.(model), data))
end

function add_poisson_noise(x::Union{Vector{Float64}, Matrix{Float64}})
  return rand.(Poisson.(x))  
end
