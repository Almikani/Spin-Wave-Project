using Sunny, GLMakie, LinearAlgebra, Random, FileIO, JLD2, Logging
using Printf: @sprintf
using Sunny: getspin
include("imports/constants.jl");
include("imports/helper_functions.jl")
include("imports/simulation.jl")

# Run with (windows)
# cd "Disordered Systems"
# $env:JULIA_NUM_THREADS=11 
# julia manual_inhom.jl

println("Threadcount: $(Threads.nthreads())")

r = range(0.1, stop=0.3, length=11)
sim_results = Array{SimulationInfo}(undef, length(r))

qpaths = [
	[[ 0,.5,     0], [ 0, 1, 0], [ 0, 2, 0]            ],
	[[ 0, 1, -0.75], [ 0, 1, 0], [ 0, 1, 1], [ 0, 1, 2]],
]

# Trying to precompile before multithreading
simulate_system(qpaths = qpaths)

Threads.@threads for i in range(1, stop=length(r))

	try
		sim_results[i] = simulate_system(
			x_ideal = r[i],
			mode = :SUN,
			units = Units(:meV, :angstrom),
			sys_repetitions = (6, 6, 6),
			qpaths = qpaths,
			print_debug = true,
			sim_id = "$(i)"
		)
	catch e
		println("Error in simulation $(i) with x_ideal=$(r[i]):\n$e")
		continue
	end

	println("Finished simulation $(i) with x_ideal=$(r[i]), x_actual=$(sim_results[i].x_actual)\n")
end


# Plotting is (apparently) not thread safe
for (i, result) in enumerate(sim_results)
	println("Saving result $(i)")

	# Inelegant i know..
	repetition_str = join(result.sys_repetitions, 'x')

	save_sim(
		result;
		save_name = string("Results/sweep SUN/$(i)_{x=$(@sprintf("%.3f", result.x_actual))}_", repetition_str, "_SUN"),
		show_plot = true,
		save_plot = true,
		save_data = true,
	)
end