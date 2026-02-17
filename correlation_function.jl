using JLD2, LaTeXStrings, Cosmology, Plots, Roots, DataStructures, Unitful, UnitfulAstro, Random, LsqFit
include("raytrace_tk.jl")


dx=1.0f0
mode="OBB"
rf=5
update=false

Hh = 0.7  # Hubble adimensional parameter
Ωm = 0.3  # Matter density parameter
Ωr = 0.0  # Radiation density parameter
ΩΛ = 0.7  # Dark energy density parameter

cosmo = Cosmology.FlatLCDM(Hh, ΩΛ, Ωm, Ωr)
Om_z( c::Cosmology.AbstractCosmology, z::Real ) = c.Ω_m * (1.0+z)^3 / Cosmology.E(c, z)^2
Om_z( c::Cosmology.FlatLCDM, z::Real ) = c.Ω_m / ( c.Ω_m + c.Ω_Λ*(1.0+z)^(-3) )
ρ_crit( c::Cosmology.AbstractCosmology, z::Real ) = 277.505425056 * (c.h*Cosmology.E(c,z))^2 #~10^9 h² Msun/Mpc³ j'ai l'impression
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

########################################### SOME NEEDED FUNCTIONS FOR FITTING ###################################################

function Flens(r, rs)
	
   	x = r / rs
	val = x^2 - 1
	
	if val > 0
		s = sqrt(val)
		return atan(s)/s
	elseif x < 1.0
		s = sqrt(abs(val))
		return atanh(s)/s
	elseif x == 1.0 
		return 1
	end
end

function Fprime(x)
	return 1/(x^3-x) * (1 - x^2*Flens(x))
end


function Hlens(x)
	return 1-x*Flens(x)-1/x
end 


function Glens(x)
	return (1-Flens(x))*(8/(x*x*x) - 20/x + 15*x)
end


function kappafunc(r, ks, rs)
	rs = rs*u"kpc"
	
	return 2 .* ks .* rs^2 .* (1 .- Flens.(r, rs)) ./ (r .^2 .- rs .^2)
end


function flexF(x)
	Fs = 1 
	return -(2*x*Flens(x)-Hlens(x))*2*Fs/(y*y-1)^2
end

function flexG(x) 
	Fs = 1
	return 2*Fs*(log(x/2)*8/(x*x*x)+((3/x)*(1-2*x*x)+Glens(x))/(x*x-1)^2)
end


#################################################################################################################################
############# END OF DEFINITION OF COSMOLOGICAL PARAMETERS AND FUNCTIONS, NOW ON TO COMPUTING CORRELATIONS ######################
#################################################################################################################################

ids = Int64[50,100,150,200,250,300,350,400,417,450]
zlist=Float64[0.21253133,0.3513419,0.5109155,0.7274435,1.0156995,1.4190572,1.8170974,2.4618596,2.79021206,3.913962228095617]

specbase = "alpha2jac_final"

arcm = 50f0
arcmdoub = 50

# 1 for jackknife, 0 for bootstrap
error_estimator = 1

# 2 for most massive galaxies, 1 for intermediate mass galaxies, 0 for low mass galaxies
massindic = 2

# 0 for central galaxies, 1 for satellite galaxies
satellite = 1



kwargs = Dict(:xlabel => "Distance", 
              :xscale => :log10, 
              :yscale => :log10, 
              :legend => true)

			  
fit_params = OrderedDict{Int64, Tuple{Vector{Float64}, Vector{Float64}}}()
m_c_dict = OrderedDict{Int64, Vector{Float64}}()

if satellite == 0
	specbase = "central_"*specbase
elseif satellite == 1
	specbase = "satellite_"*specbase
end

if massindic == 2
	specbase = "highmass_"*specbase

elseif massindic == 1
	specbase = "intermediatemass_"*specbase

elseif massindic == 0
	specbase = "lowmass_"*specbase
end



#################################################################################################################################
@assert length(ids)==length(zlist)
for i in eachindex(ids)

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
	end
	

	@assert (length(ra_img)==length(dec_img)) & (length(ra_img) == length(z)) & (length(ra_img) == length(mass))
	@info @sprintf("%d galaxies out of %d , lensing signal evaluated at z € [%.5f , %.5f], mass indicator = %d, satellite indicator = %d",num_galaxies, ngal_all, z_half-dz, z_half+dz, massindic, satellite)
	@info @sprintf("Gathered data, computing correlation function...id=%d and zs=%.5f",id,z0)
	@time Sres, rbin, kappa = comp_gs_corr(ra_img, dec_img, fn,rmax=arcm)
	@info @sprintf("Corr func computed, storing Sres and rbin in Data/corr_data/Sres_rbin_%s_id_0%d_arcm%d.jld2",specifier,id,arcmdoub)
	odfn = @sprintf("Data/corr_data/Sres_rbin_%s_id_0%d_arcm%d.jld2",specifier,id,arcmdoub)
	@save odfn Sres rbin
	GC.gc(true)
end

#################################################################################################################################
##################### CORRELATIONS COMPUTED, ON TO ESTIMATE ERRORS, FIT MODELS AND PLOT THE RESULTS #############################
#################################################################################################################################


for (i, val) in enumerate(ids)
    	z_now = zlist[i]
	comov_d = comoving_radial_dist(cosmo, z_now/2)
 
	i += 1

	specifier = specbase
	

	datafile = @sprintf("Data/corr_data/Sres_rbin_%s_id_0%d_arcm%d.jld2",specifier,val,arcmdoub)
	@load datafile Sres rbin 
	@eval $(Symbol("Sres$val")) = $Sres
	@eval $(Symbol("rbin$val")) = $rbin
	rbin_centers = 0.5*(rbin[1:end-1]+rbin[2:end])

	if get(kwargs, :yscale, :identity) == :log10
		specifier = "log_"*specifier
	else
		specifier = "nolog_"*specifier
	end

	G_temp = Sres[1,:,:] ./ Sres[5,:,:]
	F_temp = Sres[2,:,:] ./ Sres[5,:,:]
	kappa_temp = Sres[3,:,:] ./ Sres[5,:,:]
	gamma_temp = Sres[4,:,:] ./ Sres[5,:,:]

	G = mean(G_temp, dims=2)
	F = mean(F_temp, dims=2)
	gamma = mean(gamma_temp, dims=2)
	kappa = mean(kappa_temp, dims=2)
														#Check for NaN values and remove them
	nan_mask = isnan.(F) .| isnan.(G) .| isnan.(gamma) .| isnan.(kappa)
	nan_indices = findall(vec(nan_mask))	
	rbin_centers = deleteat!(rbin_centers, nan_indices)
	F = F[isfinite.(F)]
	G = G[isfinite.(G)]
	gamma = gamma[isfinite.(gamma)]
	kappa = kappa[isfinite.(kappa)]


	if error_estimator == 1 							#JACKKNIFE 
		#specifier = "jackknifeerrors_"*specifier

		n_jackknife = 5000
		n_exclude = size(Sres, 3) ÷ 100
		size1 = size(G_temp, 1)
		G_jackknife = zeros(size1, n_jackknife)
		F_jackknife = zeros(size1, n_jackknife)
		gamma_jackknife = zeros(size1, n_jackknife)
		kappa_jackknife = zeros(size1, n_jackknife)


		for i in 1:n_jackknife
			exclude_indices = Random.randperm(size(Sres, 3))[1:n_exclude]

			mask = [!(j in exclude_indices) for j in 1:size(Sres, 3)]

			G_temp_jackknife = G_temp[:, mask]
			F_temp_jackknife = F_temp[:, mask]
			gamma_temp_jackknife = gamma_temp[:, mask]
			kappa_temp_jackknife = kappa_temp[:, mask]

			G_jackknife[:, i] = mean(G_temp_jackknife, dims=2)
   			F_jackknife[:, i] = mean(F_temp_jackknife, dims=2)
    		gamma_jackknife[:, i] = mean(gamma_temp_jackknife, dims=2)
    		kappa_jackknife[:, i] = mean(kappa_temp_jackknife, dims=2)
		end
						
		mean_Gjk = mean(G_jackknife,dims=2)
		mean_Fjk = mean(F_jackknife,dims=2)
		mean_gammajk = mean(gamma_jackknife,dims=2)
		mean_kappajk = mean(kappa_jackknife,dims=2)

		G_jackknife_i = zeros(size(G_jackknife))
		F_jackknife_i = zeros(size(F_jackknife))
		gamma_jackknife_i = zeros(size(gamma_jackknife))
		kappa_jackknife_i = zeros(size(kappa_jackknife))

		for i in 1:size(Sres,3)
		    for j in 1:size(G_jackknife, 1)
				G_jackknife_i[j, i] = (G_jackknife[j, i] - mean_Gjk[j])^2
				F_jackknife_i[j, i] = (F_jackknife[j, i] - mean_Fjk[j])^2
				gamma_jackknife_i[j, i] = (gamma_jackknife[j, i] - mean_gammajk[j])^2
				kappa_jackknife_i[j, i] = (kappa_jackknife[j, i] - mean_kappajk[j])^2
		    end
		end

		G_error = sqrt.((n_jackknife - 1) * mean(G_jackknife_i, dims=2))
		F_error = sqrt.((n_jackknife - 1) * mean(F_jackknife_i, dims=2))
		gamma_error = sqrt.((n_jackknife - 1) * mean(gamma_jackknife_i, dims=2))
		kappa_error = sqrt.((n_jackknife - 1) * mean(kappa_jackknife_i, dims=2))
	
	end

	if get(kwargs, :yscale, :identity) == :log10
		@info("Absolute value for log")
		G = abs.(G)
		F = abs.(F)
		gamma = abs.(gamma)
		kappa = abs.(kappa)
	end
	
	G = abs.(G)
	F = abs.(F)
	gamma = abs.(gamma)
	kappa = abs.(kappa)

	rbin_centers = comov_d*rbin_centers*3.14159/180 
	rbin_centers = uconvert.(u"kpc", rbin_centers) 

	pG = scatter(rbin_centers, G, yerr = G_error, ylabel="G", title="G for $val"; kwargs...)
	pF = scatter(rbin_centers, F, yerr = F_error, ylabel="F", title="F for $val"; kwargs...)
	pgamma = scatter(rbin_centers, gamma, yerr = gamma_error, ylabel=L"\gamma", title=L"\gamma for"*"$val"; kwargs...)
	pkappa = scatter(rbin_centers, kappa, yerr = kappa_error, ylabel=L"\kappa", title=L"\kappa for"*"$val"; kwargs...)

	# Save the plots in the corr_plots/*lensing_component*/ directory
	allpath = @sprintf("corr_plots/all/corr_%s_%d_arcm_id_%.4d.png",specifier,arcmdoub,val)
	Gpath = @sprintf("corr_plots/G/G_%s_%d_arcm_id_%.4d.png",specifier,arcmdoub,val)
	Fpath = @sprintf("corr_plots/F/F_%s_%d_arcm_id_%.4d.png",specifier,arcmdoub,val)
	gammapath = @sprintf("corr_plots/gamma/gamma_%s_%d_arcm_id_%.4d.png",specifier,arcmdoub,val)
	kappapath = @sprintf("corr_plots/kappa/kappa_%s_%d_arcm_id_%.4d.png",specifier,arcmdoub,val)
	
	savefig(pG, Gpath)
	savefig(pF, Fpath)
	savefig(pgamma, gammapath)
	savefig(pkappa, kappapath)

	@info @sprintf("Generated correlation plots for id=%d",val)
	
	pall = scatter(rbin_centers, G, yerr = G_error, ylabel="Value", title="correlations for $val", label=L"G_+"; kwargs... )
	scatter!(rbin_centers, F, yerr = F_error, label=L"F_+"; kwargs... )
	scatter!(rbin_centers, gamma, yerr = gamma_error, label=L"\gamma"; kwargs... )
	scatter!(rbin_centers, kappa, yerr = kappa_error, label=L"\kappa"; kwargs... )

	savefig(pall,allpath)


	#Fitting model
	model(x,p) = kappafunc(x,p[1],p[2])
	p0 = [1,0.01]


	#Fitting kappa
	kappafitpath = @sprintf("corr_plots/fits/kappa/kappa_fit_%s_%d_arcm_id_%.4d.png",specifier,arcmdoub,val)
	kappaFit = curve_fit(model, rbin_centers, kappa, p0)
	kappaFitErrors = estimate_errors(kappaFit, 0.95)

	

	rvec = range(rbin_centers[1], stop=rbin_centers[end] ,length=100) 

	fittedKappa = model(rvec, collect(kappaFit.param))
	
	(mass, c) = ksrs2mc(kappaFit.param[1], kappaFit.param[2], 1.0, runit = :kpc, z=z_now/2)	
	m_c_dict[val] = [mass, c]

	fit_params[val] = (kappaFit.param, kappaFitErrors)
	if get(kwargs, :yscale, :identity) == :log10
		fittedKappa = abs.(fittedKappa)
	end

	p = scatter(rbin_centers, kappa,yerr = kappa_error, label = L"\bar\kappa", ylabel = L"\kappa", title = "Kappa fit for $val"; kwargs...)
	plot!(rvec, fittedKappa, label="NFW "*L"\kappa"*" profile fit"; kwargs...)
	savefig(p,kappafitpath)
end
