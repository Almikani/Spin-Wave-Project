# Atom types
types = ["Ni", "Fe"]

# Lattice vectors
latvecs = Dict(
    "Ni" => lattice_vectors(10.02, 5.86, 4.68, 90, 90, 90),
    "Fe" => lattice_vectors(10.337, 6.011, 4.695, 90, 90, 90),
)

# Magnetic moments
moments = Dict(
    "Ni" => Moment(s=1, g=2), 
    "Fe" => Moment(s=2, g=2),
)

# Translates a Sunny "Bond" to a letter code for the exchange interactions 
# (see Yiu et al. for details)
letter_from_bond = Dict(
    Bond(2, 3, [0, 0, 0]) => "bc",
    Bond(1, 1, [0, 0, 1]) => "c",
    Bond(2, 4, [0, 0, 0]) => "ac",
    Bond(1, 3, [0, 0, 0]) => "ac",
    Bond(3, 4, [0, 0, 0]) => "ab",
    Bond(1, 2, [0, 0, 0]) => "ab",
    Bond(1, 1, [0, 1, 0]) => "b",
)

# Exchange coupling
J = Dict(
    "NiNi" => Dict(
        "bc" =>  1.036,
        "b"  =>  0.6701,
        "c"  => -0.0469,
        "ab" =>  0.2977,
        "ac" => -0.1121,
    ),
    "FeFe" => Dict(
        "bc" => 0.77,
        "b"  => 0.30,
        "c"  => 0.14,
        "ab" => 0.14,
        "ac" => 0.05,
    ),
    "NiFe" => Dict(),  # Placeholder (will be filled in below)
)
# Adds NiFe exchange couplings to dict
J["NiFe"]["bc"] = (J["NiNi"]["bc"] + J["FeFe"]["bc"])/2
J["NiFe"]["b"]  = (J["NiNi"]["b" ] + J["FeFe"]["b" ])/2
J["NiFe"]["c"]  = (J["NiNi"]["c" ] + J["FeFe"]["c" ])/2
J["NiFe"]["ab"] = (J["NiNi"]["ab"] + J["FeFe"]["ab"])/2
J["NiFe"]["ac"] = (J["NiNi"]["ac"] + J["FeFe"]["ac"])/2

# Anisotropy
# TODO: Check if easy axis is meant to be unaligned
D = Dict(
    "Ni" => Dict(
        "a" => 0.1969,
        "b" => 0.9097,
        "c" => 0.0,
    ),
    "Fe" => Dict(
        "a" => 0.62,
        "b" => 0.00,
        "c" => 1.56,
    ),
)

