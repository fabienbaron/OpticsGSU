# Centering functions
function cog(x) #note: this takes a 2D array
    xvals = [i for i = 1:size(x, 1)]
    return vec([sum(xvals' * x) sum(x * xvals)] / sum(x))
end

function posmax(x) #note: this takes a 2D array
    pos = findmax(x)[2]
    return [pos[1], pos[2]]
end

function recenter(x, center) #note: this takes a 2D array
    Δ = [div(size(x, 1), 2) + 1, div(size(x, 2), 2) + 1] - center
    return circshift(x, round.(Δ))
end

function shift_and_add(cube)
    return reduce(+, [recenter(cube[:, :, i], posmax(cube[:, :, i])) for i = 1:size(cube, 3)])
end

function entropy(x)
    return sum(x .* log.(x .+ 1e-30))
end



#Windowing functions
function bartlett_hann1d(n, N)
    return -0.48 * abs.(n / N .- 0.5) - 0.38 * cos.(2 * pi * n / N) .+ 0.62
end

function bartlett_hann1d_centered(n, N)
    nn = n .+ (div(N + 1, 2))
    return bartlett_hann1d(nn, N)
end

function bartlett_hann2d(i, j, N)
    return bartlett_hann1d_centered(i, N) .* bartlett_hann1d_centered(j, N)
end

# function gaussian2d(n,m,sigma)
#     g2d = [exp(-((X-(m/2)).^2+(Y-n/2).^2)/(2*sigma.^2)) for X=1:m, Y=1:n]
#     return g2d
# end

function gaussian2d(n, m, sigma; normalized=false)
    g2d = [exp(-((X - (m ÷ 2 + 1)) .^ 2 + (Y - (n ÷ 2 + 1)) .^ 2) / (2 * sigma .^ 2)) for X = 1:m, Y = 1:n]
    if normalized == true
        g2d /= sum(g2d)
    end
    return g2d
end

function meshgrid(xx::Union{Vector{Float32},Vector{Float64}}) #example: meshgrid([-N/2:N/2-1;]*δ);
    x = [j for i = xx, j = xx]
    return x, x'
end

function meshgrid(xx::Int64) #example: meshgrid(N);
    x = [j for i = 1:xx, j = 1:xx] .- (div(xx, 2) + 1)
    return x, x'
end

function meshrad(xx::Union{Vector{Float32},Vector{Float64}})
    x = [j for i = xx, j = xx]
    return hypot.(x, x')
end

function meshrad(xx::Int64)
    x = [j for i = 1:xx, j = 1:xx] .- (div(xx, 2) + 1)
    return hypot.(x, x')
end

function meshpol(xx::Union{Vector{Float32},Vector{Float64}})
    x = [j for i = 1:xx, j = 1:xx]
    return hypot.(x, x'), atan.(x', x)
end

function meshpol(xx::Int64)
    x = [j for i = 1:xx, j = 1:xx] .- (div(xx, 2) + 1)
    return hypot.(x, x'), atan.(x', x)
end

function cart2pol(x, y)
    return hypot.(x, y), atan.(y, x)
end

function airmass(z) # Pickering 2002 z in rad
    h = pi / 2 - z
    return 1.0 / sin(z + 244 / (165 + 47 * h^1.1))
end
