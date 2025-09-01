"""
	lerp(x1, x2, α)

General linear interpolation between `x1` and `x2` with parameter `α ∈ [0, 1]`.
"""
@inline function lerp(x1::T, x2::T, α::Real) where T
    return (1 - α) * x1 + α * x2
end

"""
    visualize_site_types(sys::System, site_atoms::Dict)

Plots the current spin configuration, setting `Ni` to blue and `Fe` to red.
Requires `site_atoms` to be a dictionary mapping each site to its atom type.
This is not included in `Sunny` so one must be made manually like so:
```Julia
site_atoms = Dict()
for site in eachsite(sys_inhom)
	if rand() > x_ideal
		site_atoms[site] = "Ni"
	else
		site_atoms[site] = "Fe"
	end
end
```

"""
@inline function visualize_site_types(sys::System, site_atoms::Dict)
    plot_spins(sys; color=[site_atoms[site] == "Ni" ? :blue : :red for site in eachsite(sys)])
end


"""
    visualize_site_types(sys::System)

Plots the current spin configuration, setting `S < 1.5` to blue and `S > 1.5` to red.
Should correspond to blue: `Ni`, red: `Fe`.
```

"""
@inline function visualize_spin_values(sys::System)
	@inline function vector_norm(v)
		return sqrt(sum(v[i]^2 for i in eachindex(v)))
	end

    plot_spins(sys; color=[vector_norm(sys.dipoles[site]) < 1.5 ? :blue : :red for site in eachsite(sys)])
end