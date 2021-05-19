using Distributions
using Plots

population = sort(convert(Array{Int64}, rand(Bernoulli(.7), 100)))
population_over_time = []
n_iter = 100
gr()

function single_reproduction(population::AbstractArray)
    n = length(population)
    pₐ = sum(population) / length(population)
    fₐ = Bernoulli(pₐ)
    sort(convert(Array{Int64}, rand(fₐ, n)))
end

function multiple_reproductions(seed_array::AbstractArray, population_over_time::AbstractArray, n_iterations::Int64)
    seed_array = single_reproduction(seed_array)
    push!(population_over_time, seed_array)
    if n_iterations == 1
        hcat(population_over_time...)
    else
        multiple_reproductions(seed_array, population_over_time, n_iterations-1)
    end
end

function plot_many_samples(n_samples::Int64, n_iterations::Int64, Nₑ::Int64, pₐ::Float64)
    p = false
    for i in 1:n_samples
        population = sort(convert(Array{Int64}, rand(Bernoulli(pₐ), Nₑ)))
        population_at_time = multiple_reproductions(population, [], n_iterations)
        if p
            display(plot!(dropdims(sum(population_at_time, dims=1), dims=1), legend=false))
        else
            display(plot(dropdims(sum(population_at_time, dims=1), dims=1), xlims=(0, n_iterations), ylims=(0, Nₑ)))
            p = true
        end
    end
end

plot_many_samples(10, 10000, 10000, .5)
# j = multiple_reproductions(population, population_over_time, n_iter)
# plot(dropdims(sum(j, dims=1), dims=1))
# j = multiple_reproductions(population, [], n_iter)
# plot!(dropdims(sum(j, dims=1), dims=1))
# j = multiple_reproductions(population, [], n_iter)
# plot!(dropdims(sum(j, dims=1), dims=1))