using Sunny, GLMakie, LinearAlgebra, Random, FileIO, JLD2
include("imports/constants.jl");
include("imports/helper_functions.jl")
include("imports/simulation.jl")

b = load_sim_data("Disordered Systems/Results/sweep SUN/1_{x=0.000}_1x1x1_SUN")


println("Loaded simulation")
show_sim(b)