using Sunny, GLMakie, LinearAlgebra, Random
using Sunny: getspin
include("imports/constants.jl");
include("imports/helper_functions.jl")

units = Units(:meV, :angstrom);

x_ideal = 0.2  # Fe fraction
sys_repetitions = (6, 6, 6)  # Recommended (5,5,5) at most
mode = :dipole_uncorrected  # SUN:, :dipole, :dipole_uncorrected

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
x_actual = N_Fe / length(site_atoms)

println("Actual Fe fraction: $x_actual")


## REMAKE CRYSTAL WITH INTERPOLATED LATTICE VECTORS ##
lerped_latvecs = lerp(latvecs["Ni"], latvecs["Fe"], x_actual)
cryst = Crystal(lerped_latvecs, [[1/4, 1/4, 0]], 62)
sys = System(cryst, [1 => Moment(s=1, g=2)], mode)
sys_inhom = to_inhomogeneous(repeat_periodically(sys, sys_repetitions))

println("Setting interaction energies...")
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

println("Setting anisotropy energies...")
# Sets anisotropy energy and value of S for each site
for site in eachsite(sys_inhom)
	type = site_atoms[site]
	anisotropy_energy = S -> D[type]["a"]*S[1]^2 + D[type]["b"]*S[2]^2 + D[type]["c"]*S[3]^2;
	set_onsite_coupling_at!(sys_inhom, anisotropy_energy, site)
	
	# Potentially dangerous way of setting the site-specific scaling
	sys_inhom.κs[site] = moments[type].s * (sys_inhom.Ns[site]-1)/2
	set_dipole!(sys_inhom, sys_inhom.dipoles[site], site)
end

println("Minimizing energy...")
randomize_spins!(sys_inhom)
minimize_energy!(sys_inhom; maxiters=20000)


sum_Fe = [0.0, 0.0, 0.0]
sum_Ni = [0.0, 0.0, 0.0]
N_Fe = 0
N_Ni = 0

for site in eachsite(sys_inhom)
	spin = getspin(sys_inhom, site).S

	if site_atoms[site] == "Fe"
		global sum_Fe .+= spin
		global N_Fe += 1
	elseif site_atoms[site] == "Ni"
		global sum_Ni .+= spin
		global N_Ni += 1
	end
end

println("Mean Fe spin: $(sum_Fe / N_Fe)")
println("Mean Ni spin: $(sum_Ni / N_Ni)")

visualize_site_types(sys_inhom, site_atoms)