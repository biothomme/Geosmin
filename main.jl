include("./sympatricpopulation.jl")

using .SympatricPopulation


plot_population_histories(
    10, 200, [.01 .05 .1 .2 .3 .4 .5 .6 .7 .8 .9], n_chromosomes=12, 
    n_turnover=1, n_iterations=10, sample_function=crosspollination_sampling, 
    variance_plot=true, mean_heatmap=true,fixation_plot=true, 
    hybridization_content_plot=true)






