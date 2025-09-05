using Sunny, GLMakie, LinearAlgebra, Random, FileIO, JLD2, Logging

struct SimulationInfo
	# PARAMETERS
	x_ideal::Float64
	mode::Symbol
	units::Units
	sys_repetitions::Tuple{Int, Int, Int}

	kernel::Sunny.AbstractBroadening
	energies::AbstractRange
	qpaths::Vector{Sunny.QPath}

	# RESULTS
	x_actual::Float64
	N_Fe::Int
	N_Ni::Int
	crystal::Crystal
	system::System
	SWT::SpinWaveTheoryKPM
	intensities::Vector{Sunny.Intensities}
end

# function simulate_system(sim_info::SimulationInfo; print_debug::Bool = false, sim_id::Union{String, Nothing} = nothing) 
# 	return simulate_system(
# 		x_ideal = sim_info.x_ideal,
# 		mode = sim_info.mode,
# 		units = sim_info.units,
# 		sys_repetitions = sim_info.sys_repetitions,
# 		print_debug = print_debug,,
# 		sim_id = sim_id
# 	)
# end

function simulate_system(;
	x_ideal::Float64 = 0.0,
	mode::Symbol = :dipole_uncorrected,
	units::Units = Units(:meV, :angstrom),
	sys_repetitions::Tuple{Int, Int, Int} = (1, 1, 1),
	kernel::Sunny.AbstractBroadening = lorentzian(fwhm=0.4),
	energies::AbstractRange = range(0.0, 8.0, 600),
	qpaths::Vector{Vector{Vector{Float64}}},
	print_debug::Bool = false,
	sim_id::Union{String, Nothing} = nothing
)

	# Make crystal and system
	cryst = Crystal(latvecs["Ni"], [[1/4, 1/4, 0]], 62)
	sys   = System(cryst, [1 => Moment(s=1, g=2)], mode)
	sys_inhom = to_inhomogeneous(repeat_periodically(sys, sys_repetitions))

	# Populate the system randomly with Ni and Fe
	site_atoms = Dict()
	for site in eachsite(sys_inhom)
		if rand() > x_ideal
			site_atoms[site] = "Ni"
		else
			site_atoms[site] = "Fe"
		end
	end
	N_Fe = count(i->i == "Fe", values(site_atoms))
	N_Ni = count(i->i == "Ni", values(site_atoms))
	x_actual = N_Fe / (N_Fe + N_Ni)

	if print_debug println("{$sim_id} Actual Fe fraction: $x_actual") end


	## REMAKE CRYSTAL WITH INTERPOLATED LATTICE VECTORS ##
	lerped_latvecs = lerp(latvecs["Ni"], latvecs["Fe"], x_actual)
	cryst = Crystal(lerped_latvecs, [[1/4, 1/4, 0]], 62)
	sys = System(cryst, [1 => Moment(s=1, g=2)], mode)
	sys_inhom = to_inhomogeneous(repeat_periodically(sys, sys_repetitions))

	if print_debug println("{$sim_id} Setting interaction energies...") end
	for bond in keys(letter_from_bond)
		for (site1, site2, offset) in symmetry_equivalent_bonds(sys_inhom, bond)
			if site_atoms[site1] == "Ni" && site_atoms[site2] == "Ni"
				# Ni-Ni interaction
				set_exchange_at!(sys_inhom, J["NiNi"][letter_from_bond[bond]], site1, site2, offset=offset)

			elseif site_atoms[site1] == "Fe" && site_atoms[site2] == "Fe"
				# Fe-Fe interaction
				set_exchange_at!(sys_inhom, J["FeFe"][letter_from_bond[bond]], site1, site2, offset=offset)

			else
				# Ni-Fe interaction
				set_exchange_at!(sys_inhom, J["NiFe"][letter_from_bond[bond]], site1, site2, offset=offset)

			end
		end
	end

	if print_debug println("{$sim_id} Setting anisotropy energies...") end
	# Sets anisotropy energy and value of S for each site
	for site in eachsite(sys_inhom)
		type = site_atoms[site]
		anisotropy_energy = S -> D[type]["a"]*S[1]^2 + D[type]["b"]*S[2]^2 + D[type]["c"]*S[3]^2;
		set_onsite_coupling_at!(sys_inhom, anisotropy_energy, site)
		
		# Potentially dangerous way of setting the site-specific scaling
		sys_inhom.κs[site] = moments[type].s * (sys_inhom.Ns[site]-1)/2
		set_dipole!(sys_inhom, sys_inhom.dipoles[site], site)
	end

	if print_debug println("{$sim_id} Minimizing energy...") end
	randomize_spins!(sys_inhom)
	minimize_energy!(sys_inhom; maxiters=100000)

	if print_debug println("{$sim_id} Spin wave theory time...") end
	SWT = SpinWaveTheoryKPM(sys_inhom; measure=ssf_perp(sys_inhom), tol=0.05)

	# Transforms list of qpath stops into actual qpaths
	qpaths = [q_space_path(cryst, qp, 300) for qp in qpaths]

	if print_debug println("{$sim_id} Calculating intensities...") end
	intensity_list = [intensities(SWT, qpath; energies, kernel) for qpath in qpaths]

	return SimulationInfo(
		# Parameters
		x_ideal,
		mode,
		units,
		sys_repetitions,

		kernel,
		energies,
		qpaths,

		# Results
		x_actual,
		N_Fe,
		N_Ni,
		cryst,
		sys_inhom,
		SWT,
		intensity_list
	)

end

function show_sim(sim_results::SimulationInfo)
	save_sim(sim_results; save_name="", show_plot=true, save_plot=false, save_data=false)
end

function save_sim(
	sim_results::SimulationInfo;
	save_name::String,
	show_plot::Bool = false,
	save_plot::Bool = true,
	save_data::Bool = true,
	)

	# SAVING DATA
	if save_data
		try
			FileIO.mkdir(save_name)
		catch e
		end

		Base.with_logger(NullLogger()) do;
			FileIO.save(string(save_name, "/data.jld2"), "sim_info", sim_results)
		end
	end
	
	# MAKING AND SAVING FIGURE
	fig = Figure(size=(1600, 500))
	for (i, res) in enumerate(sim_results.intensities)
		ax = plot_intensities!(fig[1, i], res; sim_results.units, title="LiNi₁₋ₓFeₓPO₄ LSWT_KPM x=$(sim_results.x_actual)")
	end

	if show_plot
		display(fig, px_per_unit=1)
	end

	if save_plot
		GLMakie.save(string(save_name, "/fig.png"), fig, px_per_unit=1)
	end
end

function load_sim_data(save_name::String)
	return FileIO.load(string(save_name, "/data.jld2"), "sim_info")
end

# "Spin-Wave-Project/Disordered Systems/Results/", 