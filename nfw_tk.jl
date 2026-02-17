using Roots, Polynomials, Interpolations

δ_c( c::Real ) = c^3 / ( log(1+c)-c/(1+c) ) / 3


########################################
function mc2ksrs( mass::Real, c::Real, q::Real ; z::Real=0.3, Δ=:vir, runit=:kpc )
    (Δ==:vir) && ( Δ = Δ_vir( cosmo, z ) )
    (a2k,scrit) = dist_scaling( z )
    rho0 = ρ_crit( cosmo, z ) * Δ
    c_cor = c*sqrt(q)
    rs = ( mass /( 4π/3 * rho0 ))^(1.0/3.0) / c_cor
    ks = δ_c(c_cor) * rs * rho0/ scrit * sqrt(q)
    (runit==:arcsec) && ( rs /= a2k )
    return (ks,rs)
end
mc2ksrs( mass::Real, c::Real; z::Real=0.3, Δ=:vir, runit=:kpc ) = mc2ksrs( mass, c, 1.0; z=z, Δ=Δ, runit=runit )
function mc2ksrs( mass::Vector{Float64}, c::Vector{Float64}, q::Vector{Float64} ; z::Real=0.3, Δ=:vir, runit=:kpc )
    (Δ==:vir) && ( Δ = Δ_vir( cosmo, z ) )
    (a2k,scrit) = dist_scaling( z )
    rho0 = ρ_crit( cosmo, z ) * Δ
    c_cor = c .* sqrt.(q)
    rs = ( mass ./( 4π/3 * rho0 )).^(1.0/3.0) ./ c_cor
    ks = δ_c.(c_cor) .* rs * rho0/ scrit .* sqrt.(q)
    (runit==:arcsec) && ( rs ./= a2k )
    return (ks,rs)
end

########################################
function ksrs2mc( ks::Real, rs::Real, q::Real ; Δ=:vir, runit=:kpc, z::Real=0.3 )
    (Δ==:vir) && ( Δ = Δ_vir( cosmo, z ) )
    (a2k,scrit) = dist_scaling( z )    
    rho0 = ρ_crit( cosmo, z ) * Δ
    lrs = deepcopy(rs)
    (runit==:arcsec) && ( lrs *= a2k )
    dc = ks/sqrt(q) * scrit / rho0 / lrs
    c = find_zero( x->δ_c(x)-dc, (0.001,1000.0), xatol=1e-3 )
    mass = 4π/3 * rho0 * (c*lrs)^3
    c /= sqrt(q)
    return (mass,c)
end
ksrs2mc( ks::Real, rs::Real, Δ=:vir, runit=:kpc, z::Real=0.3 ) =
    ksrs2mc( ks, rs, 1.0, Δ=Δ, runit=runit, z=z )
#function ksrs2mc( ks::Vector{Real}, rs::Vector{Real}, q::Vector{Real} ; Δ=:vir, runit=:kpc, z::Real=0.3 )
function ksrs2mc( ks::Vector{Float64}, rs::Vector{Float64}, q::Vector{Float64} ; Δ=:vir, runit=:kpc, z::Real=0.3 )
    (Δ==:vir) && ( Δ = Δ_vir( cosmo, z ) )
    (a2k,scrit) = dist_scaling( z )    
    rho0 = ρ_crit( cosmo, z ) * Δ
    lrs = deepcopy(rs)
    (runit==:arcsec) && ( lrs .*= a2k )
    dc = ks ./ sqrt.(q) .* scrit ./ rho0 ./ lrs
    c = map( x -> find_zero( c->δ_c(c)-x, (0.001,1000.0), xatol=1e-4 ), dc )
    mass = 4π/3 * rho0 * (c.*lrs).^3
    c ./= sqrt.(q)
    return (mass,c)
end

########################################
function c_m_K16_setup()
    zt = ( [0,0.35,0.5,1.,1.44,2.15,2.50,2.90,4.10,5.40] , )
    itpc0t  = interpolate( zt, [7.4,6.25,5.65,4.3,3.53,2.70,2.42,2.20,1.92,1.65], Gridded(Linear()) )
    itpgamt = interpolate( zt, [0.120,0.117,0.115,0.110,0.095,0.085,0.080,0.080,0.080,0.080], Gridded(Linear()) )
    itplm0t = interpolate( zt, log.( [5.5e5,1.0e5,2.0e4,900,300,42,17,8.5,2.0,0.3] ), Gridded(Linear()) )
    return (itpc0t,itpgamt,itplm0t)
end

cm_itp=c_m_K16_setup()

########################################
function c_m_K16( itp, m::Real, z::Real )    ### !!! mass m must be in h^-1 Msun (m*0.7)
    x = m*1e-12
    gam = -itp[2](z)
    a = exp( -itp[3](z) )
    return  itp[1](z) * x^gam * ( 1.0 + (x*a)^0.4 )
end
c_m_K16( m::Real, z::Real ) = c_m_K16( cm_itp, m*cosmo.h, z ) ### simplified version  in Msun

########################
##nfw_priors(z::Real) = map( x->(c=poly(x);c(z)), [[-1.6024, 3.8943,-4.6164, 2.0057],
##                                                 [ 2.6126,-0.0497, 0.0693,-0.1479],
##                                                 [ 0.1205, 0.0184, 0.0368, 0.0126],
##                                                 [ 0.2166,-0.0041,-0.0555, 0.0026],
##                                                 [ 0.7201, 0.1411,-0.0186,-0.0470]] )


########################
#function nprior( ; N::Int=10_000, alpha::Real=-1.0, cmode=(:SquareU,1.0,8.0),
#                 z::Real=0.3, mrange=[1e13,5e15] , Δ=200.0 )
#
#    if abs(alpha+1.0)<1e-4
#        imr = log.(mrange)
#        m = exp.( (imr[2]-imr[1])*rand(N) .+ imr[1] )
#    else
#        imr = mrange.^(1.0+alpha)
#        m = ( (imr[2]-imr[1])*rand(N) .+ imr[1] ).^(1.0/(1.0+alpha))
#    end
#    c = (cmode[3]-cmode[2]) * rand(N) .+ cmode[2] 
#    if cmode[1]==:K16
#        itp=c_m_K16_setup()
#        c .+= map( x-> c_m_K16(itp, x, z), m )
#    end
#    (ks,rs) = mc2ksrs( m, c, ones(N); z=z, Δ=Δ, runit=:arcsec )
#    return ( m, c, ks, rs )    
#end
