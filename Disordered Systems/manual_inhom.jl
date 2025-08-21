using Sunny, GLMakie, LinearAlgebra, Random

units = Units(:meV, :angstrom);


x_ideal = 0.5         # Fe fraction
Random.seed!(0); # Rng state

latvecs = Dict()
latvecs["Ni"] = lattice_vectors(10.02, 5.86, 4.68, 90, 90, 90)
latvecs["Fe"] = lattice_vectors(10.337, 6.011, 4.695, 90, 90, 90)

types = ["Ni", "Fe"]

cryst = Crystal(latvecs["Ni"], [[1/4, 1/4, 0]], 62, types = ["Nisda"])
sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole_uncorrected)
sys_inhom = repeat_periodically(sys, (2, 2, 2))
sys_inhom = to_inhomogeneous(sys) # TODO: change to 10 10 10

#view_crystal(sys)

#randomize_spins!(sys)



Ni_J = [
    :bc =>  1.036,
    :b  =>  0.6701,
    :c  => -0.0469,
    :ab =>  0.2977,
    :ac => -0.1121,
]
Ni_D = [    
    :a => 0.1969,
    :b => 0.9097,
    :c => 0.0,
]

Ni_J = [
    :bc => 0.77,
    :b  => 0.30,
    :c  => 0.14,
    :ab => 0.14,
    :ac => 0.05,
]
Fe_D = [
    :a => 0.62,
    :b => 0.00,  # easy axis
    :c => 1.56,
]

@show sys_inhom.crystal.root
@show sys_inhom.crystal.latvecs
@show sys_inhom.crystal.recipvecs
@show sys_inhom.crystal.positions
@show sys_inhom.crystal.types
@show sys_inhom.crystal.classes
@show sys_inhom.crystal.sg
@show sys_inhom.crystal.symprec

plot_spins(sys; color=[S[3] for S in sys.dipoles])

#for site in eachsite(sys_inhom)
    #@show sys_inhom[site]
    # if sys_inhom[site] == "Ni"
    #     set_onsite_coupling_at!(site, S -> Ni_D[:a]*S[1]^2 + Ni_D[:b]*S[2]^2 + Ni_D[:c]*S[3]^2, site)
    # elseif sys_inhom[site] == "Fe"
    #     set_onsite_coupling_at!(site, S -> Fe_D[:a]*S[1]^2 + Fe_D[:b]*S[2]^2 + Fe_D[:c]*S[3]^2, site)
    # end
#end

#for bond in get_


# set_exchange!(sys, Ni_J[bc], Bond(2, 3, [0, 0, 0]))
# set_exchange!(sys, Ni_J[c ], Bond(1, 1, [0, 0, 1]))
# set_exchange!(sys, Ni_J[ac], Bond(2, 4, [0, 0, 0]))
# set_exchange!(sys, Ni_J[ac], Bond(1, 3, [0, 0, 0]))
# set_exchange!(sys, Ni_J[ab], Bond(3, 4, [0, 0, 0]))
# set_exchange!(sys, Ni_J[ab], Bond(1, 2, [0, 0, 0]))
# set_exchange!(sys, Ni_J[b ], Bond(1, 1, [0, 1, 0]))
# set_onsite_coupling!(sys, S -> Ni_D[a]*S[1]^2 + Ni_D[b]*S[2]^2 + Ni_D[c]*S[3]^2, 1)