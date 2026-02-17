include("core.jl")

using Cosmology

cosmo = cosmology( h=0.7, OmegaM=0.30, OmegaK=0.0, OmegaR=0 )
__init__() = set_cosmology( om=cosmo.Ω_m, ol=cosmo.Ω_Λ, h=cosmo.h )


############################################################
"""
    set_cosmology( ; om::Real=0.3, ol::Real=0.7 , h::Real=0.7, zl::Real=0.3, ntab::Integer=-1 )

Defines the **softlens** internal cosmological parameters
# Arguments
- `om`: omega_matter
- `ol`: omega_lambda
- `h`: Hubble parameter H_0/100 km/s/Mpc
- `zl`: lens redshift
- `ntab`: Will pre-tabulate an array of ntab values of distance ratios from the reference lens redshift to a given source redshift... It is only activated if ntab>1

"""
function set_cosmology( ; om::Real=0.3, ol::Real=0.7, h::Real=0.7, zl::Real=0.3, ntab::Integer=-1 )
    ccall( (:jsl_set_cosmology, sllib), Cvoid,
           (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int32}),
           Float64(om), Float64(ol), Float64(zl), Float64(h), Int32(ntab) )
    global cosmo = cosmology( h=h, OmegaM=om, OmegaK=(1-om-ol), OmegaR=0 )
end

############################################################
"""
    dda( z1::Real, z2::Real )

returns the angular diameter distance between two redshifts

# Arguments
- `z1`: lower redshift
- `z2`: higher redshift

"""
function dda( z1::Real, z2::Real )
    ccall( (:jsl_dda, sllib), Float64, (Ref{Float64},Ref{Float64}),
           Float64(z1), Float64(z2) )
end
############################################################
"""
    dlsdos( zl::Real, zs::Real )

returns the ratio of the angular diameter distance between lens and source to the distance between observer and source 

# Arguments
- `zl`: lens redshift
- `zs`: source redshift

"""
function dlsdos( zl::Real , zs::Real )
    ccall( (:jsl_dlsdos, sllib), Float64, (Ref{Float64}, Ref{Float64}),
           Float64(zl), Float64(zs) )
end
############################################################
"""
    dlsdos( zs::Real )

returns the ratio of the angular diameter distance between lens and source to the distance between observer and source by interpolating through a predefined table.
Lens redshift must be defined beforehand by setting cosmological parameters and tabulating the corresponding distance ratios table.

# Arguments
- `zs`: source redshift

"""
function dlsdos( zs::Real )
   ccall( (:jsl_dlsdos_t, sllib), Float64, (Ref{Float64},), Float64(zs) )
end
############################################################
"""
    dcom( z::Real )

returns the comoving distance.

# Arguments
- `z`: source redshift

"""
dcom( z::Real ) = dda( 0.0, z )*(1.0+z)
############################################################
"""
    arcsec2kpc( z::Real )

returns how many kiloparsecs one arcsec corresponds to.

# Arguments
- `z`: source redshift

"""
arcsec2kpc( z::Real ) = dda( 0.0, z )/206.26480624709637
############################################################
"""
    dlum( z::Real )

returns the luminosity distance.

# Arguments
- `z`: source redshift

"""
dlum( z::Real ) = dda( 0.0, Float64(z) )*(1.0+z)^2

########################################
"""
    dist_scaling( z::Real )

returns arcsec2kpc and critical density

# Arguments
- `z`: source redshift

"""
function dist_scaling( z::Real )
    arcsec2kpc = dda( 0.0, Float64(z) )
    sigma_crit_inf = 1.6628140e12 / arcsec2kpc
    arcsec2kpc *= 4.8481368111e-3
    return ( arcsec2kpc, sigma_crit_inf )
end

###### 

Om_z( c::Cosmology.AbstractCosmology, z::Real ) = c.Ω_m * (1.0+z)^3 / Cosmology.E(c, z)^2
Om_z( c::Cosmology.FlatLCDM, z::Real ) = c.Ω_m / ( c.Ω_m + c.Ω_Λ*(1.0+z)^(-3) )
ρ_crit( c::Cosmology.AbstractCosmology, z::Real ) = 277.505425056 * (c.h*Cosmology.E(c,z))^2

Δ_vir( c::Cosmology.AbstractCosmology, z::Real ) = ( x=Om_z(c,z)-1.0 ; 177.65287921960845 + x*(82.0-39.0*x) )

