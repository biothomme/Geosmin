module SympatricPopulation

export make_genome,
    make_diploid_genome,
    spawn_population,
    length,
    recombinate,
    reproduce,
    HybridPopulation,
    proportion_species1,
    crosspollination_sampling,
    plot_population_histories


using Distributions
using Plots
using Random
using Statistics
using StatsBase

import Base.print

# Constructors
"""
    HaploidGenome

Simulates a haploid genome as an array, each count representing a chromosome.
Each chromosome can be either `0` or `1` for two different alleles.
"""
struct HaploidGenome
    genome::Array{Int8}
end

"""
    DiploidGenome

Simulates a diploid genome as two `HaploidGenome` instances.
"""
struct DiploidGenome
    genome1::HaploidGenome
    genome2::HaploidGenome
end

"""
    HybridPopulation

Simulates a population `Nₑ` of individual `DiploidGenome` instances. 
"""
mutable struct HybridPopulation
    population::Array
end

```
    PopulationHistory

Store history of a population in a 2d array.
```
mutable struct PopulationHistory
    diary::Array{Float64,2}
end

"""
    PopulationInTime

Simulate a population over time.
Storing properties of the population (computed by `measure_function`) in 
an array of history (`PopulationHistory`)
"""
struct PopulationInTime
    population::HybridPopulation
    history::PopulationHistory
    measure_function::Function
end

"""
    VarianceAndFixationHistory

Store histories of important properties of a population.
Here we use the variance of genome assignment, loss of one or another species
and hybrid overshoot (a hybrid has at least a proportion of 0.25 of a stranger
species).
"""
mutable struct VarianceAndFixationHistory
    variances::Array{Float64, 2}
    fixations::Array{Float64, 3}
    hybrid_overshoot::Array{Float64, 2}
end

# Functions
"""
    make_genome(n_chromosomes::Int64, p_allele::Float64)

Build a `HaploidGenome` of a given number of chromosomes (`n_chromosomes`).
The assignment of each chromosome to one or the other allele is done 
by a Bernoulli trial with probability `p_allele` for `1`.
"""
function make_genome(n_chromosomes::Int64, p_allele::Float64)
    HaploidGenome(
        convert(
            Array{Int64},
            rand(Bernoulli(p_allele), n_chromosomes)
            )
        )
end

"""
    make_diploid_genome(n_chromosomes::Int64, p_allele::Float64)

Build a `DiploidGenome` with `n_chromosomes` and `p_allele` as 
2 `HaploidGenome` instances.
"""
function make_diploid_genome(n_chromosomes::Int64, p_allele::Float64)
    DiploidGenome(make_genome(n_chromosomes, p_allele),
        make_genome(n_chromosomes, p_allele))
end

"""
    spawn_population(Nₑ::Int64, n_chromosomes::Int64, p_allele::Float64)

Build a `DiploidGenome` with `n_chromosomes` and `p_allele` as 
2 `HaploidGenome` instances.
"""
function spawn_population(Nₑ::Int64, n_chromosomes::Int64, p_allele::Float64)
    population = [
        make_diploid_genome(n_chromosomes, rand(Bernoulli(p_allele), 1)[1] ? 1. : 0.) 
        for i in 1:Nₑ]
    HybridPopulation(population)
end

"""
    length(d_g::DiploidGenome)

Retrieve number of chromosomes of each `HaploidGenome` representing 
the genome length.
Both genomes should be same length.
"""
function length(d_g::DiploidGenome)
    @assert size(d_g.genome1.genome)[1] == size(d_g.genome2.genome)[1]
    size(d_g.genome1.genome)[1]
end

"""
    length(h_p::HybridPopulation)

Retrieve number of indiviuals in a `HybridPopulation` representing 
the genome length.
Both genomes should be same length.
"""
function length(h_p::HybridPopulation)
    size(h_p.population)[1]
end

Base.deepcopy(h_g::HaploidGenome) = begin
    HaploidGenome(deepcopy(h_g.genome))
end

Base.deepcopy(d_g::DiploidGenome) = begin
    DiploidGenome(deepcopy(d_g.genome1), deepcopy(d_g.genome2))
end

"""
    proportion_species1(d_g)

Retrieve proportion of genome assigned to the species encoded with 1.
"""
function proportion_species1(d_g::DiploidGenome)
    (sum(d_g.genome1.genome) + sum(d_g.genome2.genome)) / (2*length(d_g))
end


"""
    print(diploid_genome::DiploidGenome)

Print the chromosomes in given order for both genomes.
"""
function print(diploid_genome::DiploidGenome)
    println(" - Genome 1: ", diploid_genome.genome1.genome...)
    println(" - Genome 2: ", diploid_genome.genome2.genome...)
end

"""
    print(hybrid_population::HybridPopulation)

Print the genomes of specimen in a hybrid population.
"""
function print(hybrid_population::HybridPopulation)
    println("Hybrid population (",
        length(hybrid_population),
        " dipl. individuals, ", 
        length(hybrid_population.population[1]),
        " chromosomes):")
    for (i, individual) in enumerate(hybrid_population.population)
        println("- Individual ", i, ":")
        print(individual)
    end
end

"""
    recombinate(d_g::DiploidGenome, p_recombination::Float64; shuffled_array::Array=nothing)

Make random recombination events X, with X ~ Ber(`p_recombination`).
"""
function recombinate(d_g::DiploidGenome, p_recombination::Float64; shuffled_array::Array=[])
    d_g_new = deepcopy(d_g)
    if shuffled_array == []
        genome_length = length(d_g_new)
	shuffled_array = shuffle(1:genome_length)
    end
    switch = rand(Bool)
    for i in shuffled_array
        if rand(Bernoulli(p_recombination), 1)[1]
            switch = !switch
        end
        if switch
            d_g_new.genome1.genome[i], d_g_new.genome2.genome[i] = d_g_new.genome2.genome[i], d_g_new.genome1.genome[i]
        end
    end
    d_g_new
end

"""
    reproduce(parent1::DiploidGenome, parent2::DiploidGenome, p_meioticdrive::Float64)

Simulate reproduction of two diplod genomes by taking on of each parents genomes. 
`p_meioticdrive` defines the probability to choose the first genome of both parents
"""
function reproduce(parent1::DiploidGenome, parent2::DiploidGenome, p_meioticdrive::Float64)
    genome1 = rand(Bernoulli(p_meioticdrive), 1)[1] ? 
        deepcopy(parent1.genome1) : deepcopy(parent1.genome2)
    genome2 = rand(Bernoulli(p_meioticdrive), 1)[1] ? 
        deepcopy(parent2.genome1) : deepcopy(parent2.genome2)
    DiploidGenome(genome1, genome2)
end


"""
    store_history(
        population_history::PopulationHistory,
        population::HybridPopulation,
        measure_function::Function)

Store the history of a population evolving over time. 
The abstraction level of history can be defined by using `measure_function`.
"""
function store_history(
    population_history::PopulationHistory, population::HybridPopulation,
    measure_function::Function, position::Int64)
    @assert size(population_history.diary)[1] == length(population)
    new_diary_entry = sort([measure_function(d_g) for d_g in population.population])
    population_history.diary[:, position] = new_diary_entry;
end

"""
    initiate_phistory(Nₑ::Int64, n_cycles::Int64)

Make a new `PopulationHistory` instance with a given population size `Nₑ`
and `n_cycles` generations.
"""
function initiate_phistory(Nₑ::Int64, n_cycles::Int64)
    PopulationHistory(zeros(Nₑ, n_cycles))
end

"""
    population_dt(
        Nₑ::Int64, n_cycles::Int64; n_chromosomes::Int64=10, n_turnover::Int64=1,
        p_allele::Float64=0.5, p_recombination::Float64=0.5, p_meioticdrive::Float64=0.5,
        p_crosspollination::Float64=0.5, homoz_lethality::Float64=0.1, 
        measure_function::Function=proportion_species1, fitness_function::Function=genome_balance,
        sample_function::Function=random_parent_sampling)

Simulate a population over time.
History is abstracted by `measure_function` and stored.
"""
function population_dt(
    Nₑ::Int64, n_cycles::Int64; n_chromosomes::Int64=10, n_turnover::Int64=1,
    p_allele::Float64=0.5, p_recombination::Float64=0.5, p_meioticdrive::Float64=0.5,
    p_crosspollination::Float64=0.5, homoz_lethality::Float64=0.1, 
    measure_function::Function=proportion_species1, fitness_function::Function=genome_balance,
    sample_function::Function=random_parent_sampling
    )
    @assert Nₑ >= 3
    population = spawn_population(Nₑ, n_chromosomes, p_allele)
    popu_history = initiate_phistory(Nₑ, n_cycles)
    store_history(popu_history, population, measure_function, 1)
    for i in 2:n_cycles
        population = die_and_spawn(Nₑ, population, n_turnover, 
            p_recombination, p_meioticdrive, p_crosspollination,
            homoz_lethality, measure_function, fitness_function,
            sample_function)
        store_history(popu_history, population, measure_function, i)
    end
    PopulationInTime(population, popu_history, measure_function)
end

"""
    die_and_spawn(
        Nₑ::Int64, population::HybridPopulation, n_turnover::Int64, 
        p_recombination::Float64, p_meioticdrive::Float64)

Let `n_turnover` individuals of the population die and spawn new ones.
Dieing individuals are chosen randomly. Newborns are product of reproduction
of two sampled specimen (no selfing), which had recombination event
`p_recombination` as well as possible meiotic drive preferring one or the 
other genome `p_meioticdrive`.
"""
function die_and_spawn(
    Nₑ::Int64, population::HybridPopulation, n_turnover::Int64, 
    p_recombination::Float64, p_meioticdrive::Float64, p_crosspollination::Float64,
    homoz_lethality::Float64, measure_function::Function, 
    fitness_function::Function, sample_function::Function
    )
    new_population =  Array{DiploidGenome, 1}(undef, Nₑ)
    death_said = kill_the_hybrids(
        Nₑ, n_turnover, population, homoz_lethality, fitness_function)
    # chromosome_shuffle = shuffle(1:length(population.population[1]))
    for i in 1:Nₑ
        if i in death_said
            (parent1, parent2) = sample_function(population, measure_function,
                p_crosspollination)
	    parent1 = recombinate(parent1, p_recombination) # , shuffled_array=chromosome_shuffle)
	    parent2 = recombinate(parent2, p_recombination) # , shuffled_array=chromosome_shuffle)
            offspringenome = reproduce(parent1, parent2, p_meioticdrive)
            new_population[i] = offspringenome
        else
            new_population[i] = population.population[i]
        end
    end
    new_population = shuffle(new_population)
    HybridPopulation(new_population)
end

"""
    random_parent_sampling(
        population::HybridPopulation, m_f::Function, p_cp::Float64)

Sample random parent pair from population.
"""
function random_parent_sampling(
    population::HybridPopulation, m_f::Function, p_cp::Float64)
    StatsBase.sample(
                population.population, 2; replace=false)
    # if `replace == true`, one could allow selfing
end

"""
    crosspollination_sampling(
        population::HybridPopulation, measure_function::Function,
        p_crosspollination::Float64)

Sample parent pair from population like an average pollinator does.
`p_crosspollination` defines how much less the pollinators prefer the less
interesting species from the favorite species. Thus, it is always 
0 ⪓ `p_crosspollination` ≤ 1.
"""
function crosspollination_sampling(
    population::HybridPopulation, measure_function::Function,
    p_crosspollination::Float64)
    if rand(Bernoulli(.5), 1)[1]
        weights = (1).- [measure_function(d_g) for d_g in population.population]
    else
        weights = [measure_function(d_g) for d_g in population.population]
    end
    weights = [p_crosspollination + (1 - p_crosspollination)*weight for weight in weights]
    StatsBase.sample(population.population, pweights(weights), 2; replace=false)
end

"""
    kill_the_hybrids(
        Nₑ::Int64, n_turnover::Int64, population::HybridPopulation,
        homoz_lethality::Float64, measure_function::Function)

Sample individuals to die with increased hybrid lethality.
`homoz_lethality` defines the probability of killing a homozygous
non-hybrid compared to a full hybrid.
Returns array of indices of individuals to kill.
"""
function kill_the_hybrids(
    Nₑ::Int64, n_turnover::Int64, population::HybridPopulation,
    homoz_lethality::Float64, measure_function::Function)
    ε = (1 - homoz_lethality) / .5
    weights = [measure_function(individual) * ε + homoz_lethality
        for individual in population.population]
    StatsBase.sample(1:Nₑ, pweights(weights), n_turnover; replace=false)
end

"""
    genome_balance(d_g::DiploidGenome)

Compute the proportion of heterozygous genome.
The return value is min({p₁, p₀}), where p₁, p₀ are the genomic
proportion of species 1 or 0 within one individual.
"""
function genome_balance(d_g::DiploidGenome)
    f(p) = p > .5 ? 1 - p : p
    p = f(proportion_species1(d_g))
    @assert p <= .5
    @assert p >= 0.
    p
end

"""
    uniform_kill(d_g::DiploidGenome)

Shape uniform distribution to kill without hybrid disadvantage.
"""
function uniform_kill(d_g::DiploidGenome)
    1
end

"""
    repeat_population_dt(
        Nₑ::Int64, n_cycles::Int64; n_chromosomes::Int64=10, n_turnover::Int64=1,
        n_iterations::Int64=10, p_allele::Float64=0.5, p_recombination::Float64=0.5, 
        p_meioticdrive::Float64=0.5, p_crosspollination::Float64=0.5, homoz_lethality::Float64=0.1, 
        measure_function::Function=proportion_species1, fitness_function::Function=genome_balance,
        sample_function::Function=random_parent_sampling
        )

Repeat the function `population_dt` for `n_iterations` times.
Return the mean of all population histories.
"""
function repeat_population_dt(
    Nₑ::Int64, n_cycles::Int64; n_chromosomes::Int64=10, n_turnover::Int64=1,
    n_iterations::Int64=10, p_allele::Float64=0.5, p_recombination::Float64=0.5, 
    p_meioticdrive::Float64=0.5, p_crosspollination::Float64=0.5, homoz_lethality::Float64=0.1, 
    measure_function::Function=proportion_species1, fitness_function::Function=genome_balance,
    sample_function::Function=random_parent_sampling
    )
    mean_population = nothing
    properties_history = VarianceAndFixationHistory(
        Array{Float64, 2}(undef, n_iterations, n_cycles),
        Array{Float64, 3}(undef, n_iterations, 2, n_cycles),
        Array{Float64, 2}(undef, n_iterations, n_cycles))
    for i in 1:n_iterations
        population = population_dt(
            Nₑ, n_cycles; n_chromosomes=n_chromosomes, n_turnover=n_turnover,
            p_allele=p_allele, p_recombination=p_recombination, p_meioticdrive=p_meioticdrive,
            p_crosspollination=p_crosspollination, homoz_lethality=homoz_lethality, 
            measure_function=measure_function, fitness_function=fitness_function,
            sample_function=sample_function)
        if isnothing(mean_population)
            mean_population = deepcopy(population.history.diary)
        else
            mean_population = mean_population + population.history.diary
        end
        summarize_history(population.history, properties_history, i)
    end
    mean_population = (mean_population)./ n_iterations
    (mean_population, properties_history)
end

"""
    get_population_property(
        population_history::PopulationHistory, behavioural_function::Function)

Compute and return a given parameter for the current population
"""
function get_population_property(
    population_history::PopulationHistory, behavioural_function::Function)
    dropdims(mapslices(
        x -> behavioural_function(x), population_history.diary, dims=1),
        dims=1)
end

"""
    hybrid_overshoot(population_history_slice::Array)

Compute the frequency of hybrids withih a population.
"""
function hybrid_overshoot(population_history_slice::Array)
    HYBRID_CONTENT = .25  # must be < 0.5
    n_hybrids = sum((population_history_slice.>HYBRID_CONTENT).&(
            population_history_slice.<(1-HYBRID_CONTENT)))
    f_hybrids = n_hybrids / size(population_history_slice)[1]
    f_hybrids
end

"""
    lost_in_fixation(population_history_slice::Array, species1::Bool)

Compute the frequency of specimen of one species.
To be assigned to a species, an individual must share at least 0.9 of
the species genome.
"""
function lost_in_fixation(population_history_slice::Array, species1::Bool)
    HYBRID_CONTENT = .10  # must be < 0.5
    if species1
        n_homozygotes = sum((population_history_slice.>(1-HYBRID_CONTENT)))
    else
        n_homozygotes = sum((population_history_slice.<HYBRID_CONTENT))
    end
    f_homozygotes = n_homozygotes / size(population_history_slice)[1]
    f_homozygotes
end

function summarize_history(population_history::PopulationHistory,
    properties_history::VarianceAndFixationHistory, iteration_step::Int64)
    properties_history.variances[iteration_step, :] = get_population_property(
        population_history, var)
    properties_history.hybrid_overshoot[iteration_step, :] = get_population_property(
        population_history, hybrid_overshoot)
    loss_of_species = Array{Float64, 2}(undef, 2, size(population_history.diary)[2])
    for species in 0:1
        species1 = species==1
        loss_of_species[species+1, :] = get_population_property(
            population_history, x -> lost_in_fixation(x, species1))
    end
    properties_history.fixations[iteration_step, :, :] = loss_of_species
end

"""

Make different analytical plots of population histories.
"""
function plot_population_histories(
    Nₑ::Int64, n_cycles::Int64, pollination_range::Array; 
    variance_plot::Bool=false, mean_heatmap::Bool=true, fixation_plot::Bool=false,
    hybridization_content_plot::Bool=false,
    n_chromosomes::Int64=10, n_turnover::Int64=1,
    n_iterations::Int64=10, p_allele::Float64=0.5, p_recombination::Float64=0.5, 
    p_meioticdrive::Float64=0.5, homoz_lethality::Float64=0.1, 
    measure_function::Function=proportion_species1, fitness_function::Function=genome_balance,
    sample_function::Function=random_parent_sampling, postfix::String="")
    HYBRID_OVERSHOOT = .1
    gr()
    vp = variance_plot ? plot() : nothing
    fp = fixation_plot ? plot() : nothing
    hcp = hybridization_content_plot ? plot() : nothing
    mhp_array = Array{Plots.Plot, 1}(undef, size(pollination_range)[2])
    for (j, i) in enumerate(pollination_range)
        population_history, properties_history = repeat_population_dt(
            Nₑ, n_cycles; n_chromosomes=n_chromosomes, n_turnover=n_turnover, p_allele=p_allele, p_recombination=p_recombination, p_meioticdrive=p_meioticdrive, n_iterations=n_iterations,
            p_crosspollination=i, homoz_lethality=homoz_lethality, 
            measure_function=measure_function, fitness_function=fitness_function,
            sample_function=sample_function);
	variance_plot ? plot!(vp, (1:n_cycles), mapslices(sum, properties_history.hybrid_overshoot.>HYBRID_OVERSHOOT,
            dims=1)', linewidth=1, label=string(i)) : nothing
	hybridization_content_plot ? plot!(hcp, (1:n_cycles), mapslices(x -> mean(x), properties_history.variances, dims=1)',
            linewidth=.5, label=string(i)) : nothing
        fixation_plot ? plot!(fp, (1:n_cycles), dropdims(mapslices(mean,
            mapslices(x -> minimum(x), properties_history.fixations, dims=(2)),
                dims=1), dims=1)', linewidth=.5, label=string(i)) : nothing
        if mean_heatmap
            mhp_array[j] = heatmap(population_history, title=string(i), colorbar=:none)
        end
        println(string("step", i))
    end
    mhp = mean_heatmap ? plot(mhp_array..., size = (1200, 800), link = :y) : nothing
    mhp = mean_heatmap ? savefig(mhp, string("figures/heatmap_n_e", Nₑ , postfix, ".svg")) : nothing 
    fixation_plot ? savefig(
        fp, string("figures/speciesextinction_n_e", Nₑ , postfix, ".svg")) : nothing
    variance_plot ? savefig(
        vp, string("figures/variance_n_e", Nₑ , postfix, ".svg")) : nothing
    hybridization_content_plot ? savefig(
        hcp, string("figures/hybrid_fraction_n_e", Nₑ , postfix, ".svg")) : nothing
end
# Main



end  # module
