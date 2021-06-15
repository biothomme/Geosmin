include("./sympatricpopulation.jl")

using .SympatricPopulation


plot_population_histories(
    1000, 2000, [.01 .05 .1 .2 .3 .4 .5 .6 .7 .8 .9 1. .21 .67], n_chromosomes=62÷2, 
    n_turnover=10, n_iterations=100, sample_function=crosspollination_sampling, 
    variance_plot=true, mean_heatmap=true,fixation_plot=true, 
    hybridization_content_plot=true, postfix="primula")
plot_population_histories(
    1000, 2000, [.01 .05 .1 .2 .3 .4 .5 .6 .7 .8 .9 1. .30 .42], n_chromosomes=26÷2,
    n_turnover=10, n_iterations=100, sample_function=crosspollination_sampling,
    variance_plot=true, mean_heatmap=true,fixation_plot=true,
    hybridization_content_plot=true, postfix="rhododendron")
plot_population_histories(
    1000, 2000, [.01 .05 .1 .2 .3 .4 .5 .6 .7 .8 .9 1. .39 .85 .75], n_chromosomes=36÷2,
    n_turnover=10, n_iterations=100, sample_function=crosspollination_sampling,
    variance_plot=true, mean_heatmap=true,fixation_plot=true,
    hybridization_content_plot=true, postfix="gentiana")












