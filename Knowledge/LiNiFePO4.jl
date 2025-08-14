using Sunny, GLMakie, Random, DelimitedFiles
units = Units(:meV, :angstrom);
latvecs = lattice_vectors(10.02, 5.86, 4.68, 90, 90, 90);
cryst = Crystal(latvecs, [[1/4, 1/4, 0]], 62, types = ["Ni"]);
sys = System(cryst, [1 => Moment(s=1, g=2)], :dipole_uncorrected;); # initiate all sites with Ni with S = 1, use :dipole_uncorrected for S -> infinity limit
Jbc = 1;
Jab = 0.3*Jbc;
DNi = [0.3 1.8 0];
DFe = [0.6 0 1.6];
D = [DNi; DFe];
#set_exchange!(sys, Jbc, Bond(3, 2, [0,0,0]))
#set_exchange!(sys, Jab, Bond(1, 2, [0,0,0]))
#set_exchange!(sys, Jab, Bond(3, 4, [0,0,0]))
#set_onsite_coupling!(sys, S -> DNi[1]*S[1]^2 + DNi[2]*S[2]^2 + DNi[3]*S[3]^2, 1)
#randomize_spins!(sys)
#minimize_energy!(sys)
#view_crystal(sys)
#swt = SpinWaveTheory(sys; measure=ssf_perp(sys))
#res = intensities(swt, path; energies, kernel)
#plot_intensities(res)

# large system 
L = 10; N = 4*L*L*L;
x = 0.2;
sys_inhom = to_inhomogeneous(repeat_periodically(sys, (L,L,L)))
randomize_spins!(sys_inhom)

# assigning randomly Ni or Fe to all sites and add the corresponding anitropy
flavor = Array{Int64}(undef, L,L,L,4);
for i in CartesianIndices(flavor)
    tmp = rand();
    if tmp <= x
        flavor[i] = 2; # Fe
    else
        flavor[i] = 1; # Ni
    end
    set_onsite_coupling_at!(sys_inhom, S -> flavor[i]*flavor[i]*(D[flavor[i],1]*S[1]^2 + D[flavor[i],2]*S[2]^2 + D[flavor[i],3]*S[3]^2), i)
end

# loop over bonds and assign the correct exchange (factor 1, 2 and 4 for Ni-Ni, Ni-Fe and Fe-Fe bonds respectively)
for (site1, site2, offset) in symmetry_equivalent_bonds(sys_inhom, Bond(3,2,[0,0,0]))
    set_exchange_at!(sys_inhom, flavor[site1]*flavor[site2]*Jbc, site1, site2; offset)
end
for (site1, site2, offset) in symmetry_equivalent_bonds(sys_inhom, Bond(1,2,[0,0,0]))
    set_exchange_at!(sys_inhom, flavor[site1]*flavor[site2]*Jab, site1, site2; offset)
end
for (site1, site2, offset) in symmetry_equivalent_bonds(sys_inhom, Bond(3,4,[0,0,0]))
    set_exchange_at!(sys_inhom, flavor[site1]*flavor[site2]*Jab, site1, site2; offset)
end

# find ground state
minimize_energy!(sys_inhom)

#plot_spins(sys_inhom; color=[S[3] for S in sys_inhom.dipoles], show_cell=false, ndims=3)
plot_spins(sys_inhom; color=flavor, show_cell=false, ndims=3)

# calculate spectrum
#qs = [[0, 0.5, 0], [0, 1.5, 0]]; # along (0,K,0)
qs = [[0, 0.5, -0.44], [0, 1.5, -0.44]]; # along (0,K,0)
#qs = [[0, 1, -0.5], [0, 1, 0.5]]; # along (0,1,L)
path = q_space_path(cryst, qs, 51;)
kernel = lorentzian(fwhm=0.5);
energies = range(0.0, 8.0, 50);
swt = SpinWaveTheoryKPM(sys_inhom; measure=ssf_perp(sys_inhom), tol=0.01);
res = intensities(swt, path; energies, kernel);
plot_intensities(res; units)

#writedlm("/home/elfogh/Documents/DTUprojects/LithiumOrthophosphates/LiNiFePO4/sunny/LiNi0.8Fe0.2PO4_Q=(0,qk,-0.44)_L=10_1.txt", res.data)
#writedlm("/home/elfogh/Documents/DTUprojects/LithiumOrthophosphates/LiNiFePO4/sunny/LiNi0.8Fe0.2PO4_Q=(0,K,0)_L=10_3.txt", res.data)
#writedlm("/home/elfogh/Documents/DTUprojects/LithiumOrthophosphates/LiNiFePO4/sunny/LiNi1.0Fe0.0PO4_Q=(0,K,0)_L=10_2.txt", res.data)


#how to access the individual bonds and sites?
#use length() to figure out how far to loop?
#it looks like each site is actually specified by 4 indices: x,y,z specifying which unit cell and then n specifying which ion inside the cell


#include also Jb eventually, and maybe the rest of them for completeness although it will quickly get complicated


