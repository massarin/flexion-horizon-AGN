using JLD2, LaTeXStrings, Cosmology, Plots, Random, LsqFit, Unitful, UnitfulAstro
include("raytrace_tk.jl")



dx=1.0f0
mode="OBB"
rf=5
update=false

H0 = 70  # Hubble constant in km/s/Mpc
Ωm = 0.3  # Matter density parameter
Ωr = 0.0  # Radiation density parameter
ΩΛ = 0.7  # Dark energy density parameter

cosmo = Cosmology.FlatLCDM(H0, Ωm, Ωr, ΩΛ)

Om_z( c::Cosmology.AbstractCosmology, z::Real ) = c.Ω_m * (1.0+z)^3 / Cosmology.E(c, z)^2
Om_z( c::Cosmology.FlatLCDM, z::Real ) = c.Ω_m / ( c.Ω_m + c.Ω_Λ*(1.0+z)^(-3) )
ρ_crit( c::Cosmology.AbstractCosmology, z::Real ) = 277.505425056 * (c.h*Cosmology.E(c,z))^2
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
    c = find_zero( x->δ_c(x)-ustrip(dc), (0.001,1000.0), xatol=1e-3 )
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
	return 1-x*f(x)-1/x
end 


function Glens(x)
	return (1-f(x))*(8/(x*x*x) - 20/x + 15*x)
end
####################################################

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

####################################################


####################################################




#1 for Jackknife 
#0 for Bootstrap
error_estimator = 1


kwargs = Dict(:xlabel => "Distance", 
              :xscale => :log10, 
              :yscale => :log10, 
              :legend => true)

			  
fit_params = Dict{Int64, Tuple{Vector{Float64}, Vector{Float64}}}()
m_c_dict = Dict{Int64, Vector{Float64}}()
arcm = 50f0
arcmdoub = 50

for (i, val) in enumerate(ids)
    z_now = zlist[i]
	comov_d = comoving_radial_dist(cosmo, z_now/2)
 
	i += 1

	specifier = "alpha2jacwcrossterms"
	

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
		specifier = "jackknifeerrors_"*specifier

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
	
	pall = plot(rbin_centers, G, yerr = G_error, ylabel="Value", title="correlations for $val", label=L"G_+"; kwargs... )
	plot!(rbin_centers, F, yerr = F_error, label=L"F_+"; kwargs... )
	plot!(rbin_centers, gamma, yerr = gamma_error, label=L"\gamma"; kwargs... )
	plot!(rbin_centers, kappa, yerr = kappa_error, label=L"\kappa"; kwargs... )

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
println(fit_params)


cvec = range(0.1, stop=1.0, length=100)
"""
To keep for future reference (bootstrap tentative implementation)

elseif error_estimator == 0							#BOOTSTRAP
	specifier = "bootstraperrors_"*specifier
	n_bootstrap = size(Sres, 3)
	G_bootstrap = F_bootstrap = gamma_bootstrap = kappa_bootstrap = []
	#BOOTSTRAP
	for i in 1:n_bootstrap

			  bootstrap_sample = rand(1:size(Sres, 3), size(Sres, 3))

			G_temp_bootstrap = mean(Sres[1, :, bootstrap_sample] ./ Sres[5, :, bootstrap_sample], dims=2)
		F_temp_bootstrap = mean(Sres[2, :, bootstrap_sample] ./ Sres[5, :, bootstrap_sample], dims=2)
		 gamma_temp_bootstrap = mean(Sres[4, :, bootstrap_sample] ./ Sres[5, :, bootstrap_sample], dims=2)
		kappa_temp_bootstrap = mean(Sres[3, :, bootstrap_sample] ./ Sres[5, :, bootstrap_sample], dims=2)

		push!(G_bootstrap, G_temp_bootstrap[1])
		push!(F_bootstrap, F_temp_bootstrap[1])
		push!(gamma_bootstrap, gamma_temp_bootstrap[1])
		push!(kappa_bootstrap, kappa_temp_bootstrap[1])
	end

	# Calculate bootstrap errors
	G_error = sqrt(var(G_bootstrap))
	F_error = sqrt(var(F_bootstrap))
	gamma_error = sqrt(var(gamma_bootstrap))
	kappa_error = sqrt(var(kappa_bootstrap))	


"""