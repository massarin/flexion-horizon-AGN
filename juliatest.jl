using JLD2, LaTeXStrings, Cosmology, Plots, Random, LsqFit, Unitful, UnitfulAstro, DataStructures, Roots
include("raytrace_tk.jl")

### HERE JUST THE SECOND PART OF correlation_function.jl, if you have the correlations computed and only need some plots 

dx=1.0f0
mode="OBB"
rf=5
update=false

Hh = 0.7  # Hubble parameter, adimensional
Ωm = 0.3  # Matter density parameter
Ωr = 0.0  # Radiation density parameter
ΩΛ = 0.7  # Dark energy density parameter

cosmo = Cosmology.FlatLCDM(Hh, ΩΛ, Ωm, Ωr)

Om_z( c::Cosmology.AbstractCosmology, z::Real ) = c.Ω_m * (1.0+z)^3 / Cosmology.E(c, z)^2
Om_z( c::Cosmology.FlatLCDM, z::Real ) = c.Ω_m / ( c.Ω_m + c.Ω_Λ*(1.0+z)^(-3) )
ρ_crit( c::Cosmology.AbstractCosmology, z::Real ) = 277.505425056 * (c.h*Cosmology.E(c,z))^2 #~10^9 h² Msun/Mpc³ i think ?? not my code so..
δ_c( c::Real ) = c^3 / ( log(1+c)-c/(1+c) ) / 3

Δ_vir( c::Cosmology.AbstractCosmology, z::Real ) = ( x=Om_z(c,z)-1.0 ; 177.65287921960845 + x*(82.0-39.0*x) )

function dist_scaling( z::Real )
	#arcsec2kpc = dda( 0.0, Float64(z) )
	arcsec2kpc =  angular_diameter_dist(cosmo, z)
	arcsec2kpc = uconvert(u"kpc", arcsec2kpc)
	sigma_crit_inf = 1.6628140e12 / arcsec2kpc
	#arcsec2kpc *= 4.8481368111e-3
    return ( arcsec2kpc, sigma_crit_inf )
end

function ksrs2mc( ks::Real, rs::Real, q::Real ; Δ=:vir, runit=:kpc, z::Real=0.3 )
	(Δ==:vir) && ( Δ = Δ_vir( cosmo, z ) )
	(a2k,scrit) = dist_scaling( z )    
	rho0 = ρ_crit( cosmo, z ) * Δ
	lrs = deepcopy(rs)
	(runit==:arcsec) && ( lrs *= a2k )
	dc = ks/sqrt(q) * scrit / rho0 / lrs
	c = find_zero( x->δ_c(x)-ustrip(dc), (0.001,10000.0), xatol=1e-3 )
	mass = 4π/3 * rho0 * (c*lrs)^3
	c /= sqrt(q)
    return (mass,c)
end

ids = Int64[50,100,150,200,250,300,350,400,417,450]
zlist=Float64[0.21253133,0.3513419,0.5109155,0.7274435,1.0156995,1.4190572,1.8170974,2.4618596,2.79021206,3.913962228095617]


#ids=reverse(ids)
#zlist=reverse(zlist)



####################################################
function atanhf(x)
	return 1-2*atanh(sqrt((1-x)/(1+x)))/sqrt(1-x*x)
end


function atanf(x)
	return 1-2*atan(sqrt((x-1)/(1+x)))/sqrt(x*x-1)
end

function Flens(r, rs)
	
   	x = r / rs
	val = x^2 - 1

	if val > 0
		s = sqrt(x^2 - 1)
		return 1-2*atan(sqrt((x-1)/(x+1)))/s
	elseif val < 0
		s = sqrt(-x^2 + 1)
		return 1-2*atanh(sqrt((1-x)/(1+x)))/s
	elseif val == 0
		return 1
	end
end

function Fprime(r, rs)
	x = r / rs
	return 1/(x^3-x) * (1 - x^2*Flens(r,rs))
end


function Hlens(r,rs)
	y = r/rs	
	
	return 1-y*Flens(r,rs)-1/y
end 


function Glens(r,rs)
	y = r/rs
	return (1-Flens(r,rs))*(8/(y*y*y) - 20/y + 15*y)
end

function gammaLens(r,rs)
	y = r/rs
	check = y-1
	if check < 0
		return 8*atanh(sqrt((1-y)/(1+y)))/(y*y*sqrt(1-y*y)) + 4*log(y/2)/(y*y) - 2/(y*y-1) + 4*atanh(sqrt((1-y)/(1+y)))/((y*y-1)*sqrt(1-y*y))
	elseif check > 0
		return 8*atan(sqrt((y-1)/(y+1)))/(y*y*sqrt(y*y-1)) + 4*log(y/2)/(y*y) - 2/(y*y-1) + 4*atan(sqrt((y-1)/(y+1)))/((y*y-1)*sqrt(y*y-1))
	elseif check == 0
		return 10/3 + 4*log(1/2)
	end
end
####################################################

function kappafunc_nfw(r, rs, ks)
	rs = rs*u"kpc"
	
	return 2 .* ks .* rs^2 .* (Flens.(r, rs)) ./ (r .^2 .- rs .^2)
end

function gammafunc_nfw(r,rs,gs)
	rs = rs*u"kpc"
	invGs = (1/gs)*u"kpc"
	return rs ./ invGs .* gammaLens.(r,rs)
end

function flexF_nfw(r,rs,fs)
	rs = rs*u"kpc"
	y = r/rs
	return -(2 .* y .* Flens.(r,rs) .- Hlens.(r,rs)).* 2*fs ./(y .^2 .-1).^2
end

function flexG_nfw(r,rs,fs) 
	rs = rs*u"kpc"
	y = r/rs
	return 2*fs .* (( log.(y/2) .* 8 ./ (y .* y .* y)) .+ ((3 ./ y) .* (1 .- 2 .* y .* y) .+ Glens.(r,rs)) ./ ((y .* y .- 1) .^2))
end

####################################################


####################################################



specbase = "alpha2jac_corrected"

arcm = 50f0
arcmdoub = 50

# 1 for jackknife, 0 for bootstrap
error_estimator = 1

# 2 for most massive galaxies, 1 for intermediate mass galaxies, 0 for low mass galaxies
massindic = 0

# 0 for central galaxies, 1 for satellite galaxies
satellite = 0



kwargs = Dict(:xlabel => "Distance", 
              :xscale => :log10, 
              :yscale => :identity, 
              :legend => true)

			  
kappa_fit_params = OrderedDict{Int64, Tuple{Vector{Float64}, Vector{Float64}}}()
gamma_fit_params = OrderedDict{Int64, Tuple{Vector{Float64}, Vector{Float64}}}()
F_fit_params = OrderedDict{Int64, Tuple{Vector{Float64}, Vector{Float64}}}()
G_fit_params = OrderedDict{Int64, Tuple{Vector{Float64}, Vector{Float64}}}()
m_c_dict = OrderedDict{Int64, Vector{Float64}}()


most_massive_mean_mass = OrderedDict()
least_massive_mean_mass = OrderedDict()

if satellite == 0
	specbase = "central_"*specbase
	galtitle = "central"
elseif satellite == 1
	specbase = "satellite_"*specbase
	galtitle = "satellite"
end

if massindic == 2
	specbase = "highmass_"*specbase
	masstitle = "most massive"

elseif massindic == 1
	specbase = "intermediatemass_"*specbase
	masstitle = "intermediate mass"

elseif massindic == 0
	specbase = "lowmass_"*specbase
	masstitle = "least massive"
end



@info @sprintf("Starting analysis for the %s %s galaxies",masstitle,galtitle)
for (i, val) in enumerate(ids)

	specifier = specbase	

	z0 = zlist[i]
	dz = 0.1
	id = ids[i]

	half_d = 0.5*comoving_radial_dist(cosmo, z0)
	
	f(z) = comoving_radial_dist(cosmo,z) - half_d
	
	z_half = find_zero(f,z0/2)


	fn = full_bin_name( id, mode=mode, size=dx )
	#Getting useful data from catalog
	fitsfile = FITS("Data/Galaxies_0-6_lensed.v2.0_cut_i27.fits")
	z = read(fitsfile[2],"z_true")
	ra_img = read(fitsfile[2],"RA_IMG")
	dec_img = read(fitsfile[2],"DEC_IMG")
	mass = read(fitsfile[2],"MTOTH_MAIN")
	masssubhalo = read(fitsfile[2],"MTOTH_SUB")
	halo_id = read(fitsfile[2],"IDH_MAIN")
	subhalo_id = read(fitsfile[2],"IDH_SUB")
	close(fitsfile)

	#Selecting by redshift
	mask = (z .>= (z_half - dz)) .& (z .<= (z_half + dz))
	ra_img = ra_img[mask]
	dec_img = dec_img[mask]
	mass = mass[mask]
	masssubhalo = masssubhalo[mask]
	halo_id = halo_id[mask]
	subhalo_id = subhalo_id[mask]	
	z = z[mask]

	#Satellite or central galaxies
	if satellite == 1
		sat_mask = (subhalo_id .!= halo_id)
		ra_img = ra_img[sat_mask]
		dec_img = dec_img[sat_mask]
		mass = mass[sat_mask]
		masssubhalo = masssubhalo[sat_mask]
		z = z[sat_mask]
	elseif satellite == 0
		central_mask = (subhalo_id .== halo_id)
		ra_img = ra_img[central_mask]
		dec_img = dec_img[central_mask]
		mass = mass[central_mask]
		masssubhalo = masssubhalo[central_mask]
		z = z[central_mask]
	end
	low_quantile = 0.4
	high_quantile = 0.6
	quantile_dif = 0.001
	low_mass_threshold = quantile(mass, low_quantile) 
	high_mass_threshold = quantile(mass, high_quantile) 

	#Highest masses/Intermediate/Lowest masses
	if massindic == 2
		descending_mass_indices = sortperm(mass, rev=true)
		ngal_all = length(descending_mass_indices)
		num_galaxies = min(500, ngal_all)
		massive_galaxies = descending_mass_indices[1:num_galaxies]
		ra_img = ra_img[massive_galaxies]
		dec_img = dec_img[massive_galaxies]
		z = z[massive_galaxies]
		mass = mass[massive_galaxies]
		most_massive_mean_mass[val] = mean(mass)
	elseif massindic == 1
		median_mass = median(mass)
		diff_from_median = abs.(mass .- median_mass)
		intermediate_masses_ind = sortperm(diff_from_median)
		ngal_all = length(intermediate_masses_ind)
		num_galaxies = min(500, ngal_all)
		intermediate_masses = intermediate_masses_ind[1:num_galaxies]
		ra_img = ra_img[intermediate_masses]
		dec_img = dec_img[intermediate_masses]
		z = z[intermediate_masses]
		mass = mass[intermediate_masses]
	elseif massindic == 0
		ascending_mass_indices = sortperm(mass)
		ngal_all = length(ascending_mass_indices)
		num_galaxies = min(500, ngal_all)
		least_massive_galaxies = ascending_mass_indices[1:num_galaxies]
		ra_img = ra_img[least_massive_galaxies]
		dec_img = dec_img[least_massive_galaxies]
		z = z[least_massive_galaxies]
		mass = mass[least_massive_galaxies]
		least_massive_mean_mass[val] = mean(mass)
	end
	
end
println(most_massive_mean_mass)
println(least_massive_mean_mass)

