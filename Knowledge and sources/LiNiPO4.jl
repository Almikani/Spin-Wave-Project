using Sunny, GLMakie
units = Units(:meV, :angstrom);
latvecs = lattice_vectors(10.02, 5.86, 4.68, 90, 90, 90);
positions = [[1/4,1/4,0]];
cryst = Crystal(latvecs, positions, 62)
#view_crystal(cryst)
sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole)
Jbc = 1.04;
Jb  = 0.670;
Jc  = -0.05;
Jab = 0.30;
Jac = -0.11;
Da  = 0.339;
Db  = 1.82;
set_exchange!(sys, Jbc, Bond(2, 3, [0, 0, 0]))
set_exchange!(sys, Jc, Bond(1, 1, [0, 0, -1])) # in a different unit cell
set_exchange!(sys, Jb, Bond(1, 1, [0, 1, 0]))
set_exchange!(sys, Jab, Bond(1, 2, [0, 0, 0]))
set_exchange!(sys, Jab, Bond(3, 4, [0, 0, 0]))
set_exchange!(sys, Jac, Bond(3, 1, [0, 0, 0]))
set_exchange!(sys, Jac, Bond(4, 2, [0, 0, 0]))
set_onsite_coupling!(sys, S -> Da*S[1]^2 + Db*S[2]^2, 1)
view_crystal(sys)
randomize_spins!(sys)
minimize_energy!(sys)
plot_spins(sys; color=[S[3] for S in sys.dipoles])
swt = SpinWaveTheory(sys; measure=ssf_perp(sys))
qs = [[0, 1, 0], [2, 1, 0]]
path = q_space_path(cryst, qs, 400)
res = intensities_bands(swt, path)
plot_intensities(res; units, title="LiNiPO4 LSWT")