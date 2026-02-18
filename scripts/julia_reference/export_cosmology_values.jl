# Export cosmology values from Julia for validation against Python

include("../../cosmo.jl")  # Load Julia cosmology functions
include("../../core.jl")    # Load core functions with jsl_dda

using Printf, DelimitedFiles

# Test redshifts
z_test = [0.0, 0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0]

# Output file
outfile = "scripts/julia_reference/cosmology_values.csv"

# Open file and write header
open(outfile, "w") do io
    println(io, "z_lens,z_source,D_A_lens_Mpc,D_A_source_Mpc,D_A_ls_Mpc,sigma_crit_Msun_pc2")
    
    # Loop over lens and source redshifts
    for z_lens in z_test
        for z_source in z_test
            if z_source <= z_lens
                continue  # Source must be behind lens
            end
            
            # Angular diameter distances using Julia cosmology
            # jsl_dda(z) returns D_A in Mpc
            D_A_lens = jsl_dda(z_lens)
            D_A_source = jsl_dda(z_source)
            D_A_ls = jsl_dda(z_source, z_lens)  # From lens to source
            
            # Critical surface density
            # Σ_crit = (c²/4πG) * (D_s / D_l D_ls)
            # In units of Msun/pc²
            c_km_s = 299792.458  # Speed of light km/s
            G_pc_Msun = 4.30091e-3  # Gravitational constant in pc (km/s)² Msun⁻¹
            
            if D_A_ls > 0
                sigma_crit = (c_km_s^2 / (4π * G_pc_Msun)) * (D_A_source / (D_A_lens * D_A_ls)) / 1e12  # Mpc² -> pc²
            else
                sigma_crit = Inf
            end
            
            @printf(io, "%.6f,%.6f,%.6f,%.6f,%.6f,%.6e\n",
                   z_lens, z_source, D_A_lens, D_A_source, D_A_ls, sigma_crit)
        end
    end
end

println("Cosmology values exported to $outfile")
println("Test with: julia scripts/julia_reference/export_cosmology_values.jl")
