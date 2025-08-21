using Sunny, GLMakie

units = Units(:meV, :angstrom);

latvecs = Dict()
latvecs["Ni"] = lattice_vectors(10.02, 5.86, 4.68, 90, 90, 90)
latvecs["Fe"] = lattice_vectors(10.337, 6.011, 4.695, 90, 90, 90)

println("Creating Crystal...")

cryst = Crystal(latvecs["Ni"], [[1/4, 1/4, 0]], 62, types = ["Ni"])
sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole_uncorrected)

println("Settin up exchange interactions...")

Jbc =  1.036
Jb  =  0.6701
Jc  = -0.0469
Jac = -0.1121
Jab =  0.2977
Da  =  0.1969
Db  =  0.9097
Dc  =  0.0
set_exchange!(sys, Jbc, Bond(2, 3, [0, 0, 0]))
set_exchange!(sys, Jc , Bond(1, 1, [0, 0,-1]))
set_exchange!(sys, Jb , Bond(1, 1, [0, 1, 0]))
set_exchange!(sys, Jab, Bond(1, 2, [0, 0, 0]))
set_exchange!(sys, Jab, Bond(3, 4, [0, 0, 0]))
set_exchange!(sys, Jac, Bond(3, 1, [0, 0, 0]))
set_exchange!(sys, Jac, Bond(4, 2, [0, 0, 0]))
set_onsite_coupling!(sys, S -> Da*S[1]^2 + Db*S[2]^2 + Dc*S[3]^2, 1)

println("Minimizing energy...")

randomize_spins!(sys)
minimize_energy!(sys)
plot_spins(sys; color=[S[3] for S in sys.dipoles])

println("Setting up qpath...")

qs = [[0, 0.5, 0], [0, 2, 0]]
path = q_space_path(cryst, qs, 600)
println("Performing spin wave calculations...")

kernel = lorentzian(fwhm=0.4)
energies = range(0.0, 12.0, 150)
swt = SpinWaveTheory(sys; measure=ssf_perp(sys))
res = intensities(swt, path; energies, kernel)
plot_intensities(res)

# Part 2
println("Setting up inhomogeneous system...")

sys_inhom = to_inhomogeneous(repeat_periodically(sys, (2, 2, 2))) # TODO: change to 10 10 10

for (index, (site1, site2, offset)) in enumerate(symmetry_equivalent_bonds(sys_inhom, Bond(1, 1, [1, 0, 0])))
	println(index, site1, site2, offset)
	noise = 0.0
	set_exchange_at!(sys_inhom, 1.0 + noise, site1, site2; offset)
end

println("Minimizing energy of inhomogeneous system...")

randomize_spins!(sys_inhom)
minimize_energy!(sys_inhom, maxiters=10_000)
plot_spins(sys_inhom; color=[S[3] for S in sys_inhom.dipoles], ndims=3)

println("Performing spin wave calculations on inhomogeneous system...")

swt = SpinWaveTheoryKPM(sys_inhom; measure=ssf_perp(sys_inhom), tol=0.05)
res = intensities(swt, path; energies, kernel)
plot_intensities(res)