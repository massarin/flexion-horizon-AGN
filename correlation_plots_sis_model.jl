using JLD2, LaTeXStrings, Cosmology, Plots, Random, LsqFit, Unitful, UnitfulAstro, DataStructures, Roots
include("raytrace_tk.jl")
include("nfw_tk.jl")

### HERE JUST THE SECOND PART OF correlation_function.jl, if you have the correlations computed and only need some plots 

dx=1.0f0
mode="OBB"
rf=5
update=false

Hh = 0.704  # Hubble parameter, adimensional
Ωm = 0.272  # Matter density parameter
Ωr = 0.0  # Radiation density parameter
ΩΛ = 0.728  # Dark energy density parameter

cosmo = cosmology(;h = Hh, OmegaK = 0, OmegaM = Ωm, OmegaR = nothing)

Om_z( c::Cosmology.AbstractCosmology, z::Real ) = c.Ω_m * (1.0+z)^3 / Cosmology.E(c, z)^2
Om_z( c::Cosmology.FlatLCDM, z::Real ) = c.Ω_m / ( c.Ω_m + c.Ω_Λ*(1.0+z)^(-3) )
ρ_crit( c::Cosmology.AbstractCosmology, z::Real ) = 277.505425056 * (c.h*Cosmology.E(c,z))^2 #~10^9 h² Msun/Mpc³ i think ?? not my code so..
δ_c( c::Real ) = (1/3) * c^3 / ( log(1+c)-c/(1+c) )

Δ_vir( c::Cosmology.AbstractCosmology, z::Real ) = ( x=Om_z(c,z)-1.0 ; 177.65287921960845 + x*(82.0-39.0*x) )

function dist_scaling( z::Real )
	#arcsec2kpc = dda( 0.0, Float64(z) )
	arcsec2kpc =  angular_diameter_dist(cosmo, z)
	arcsec2kpc = uconvert(u"pc", arcsec2kpc)
	arcsec2kpc = ustrip(arcsec2kpc)
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
	c = find_zero( x->δ_c(x)-dc, (1e-3,1e10), xatol=1e-3 )
	mass = 4π/3 * rho0 * (c*lrs)^3
	c /= sqrt(q)
    return (mass,c)
end


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

cm_itp=c_m_K16_setup()
function m2ksrs(cm_itp, m , z)
	q = 1
	c= c_m_K16(cm_itp, m, z)
	return mc2ksrs(m, c, q; z)
end




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
	val = y^2 - 1
	val2 = y - 1
	if val > 0 	&& val2 > 0
		s = sqrt(y^2 - 1)
		return 2*y*atan(sqrt((y-1)/(y+1)))/s - 1/y
	elseif val < 0 && val2 < 0
		s = sqrt(-y^2 + 1)
		return 2*y*atanh(sqrt((1-y)/(1+y)))/s - 1/y
	elseif val == 0
		return 0
	end

end 


function Glens(r,rs)
	y = r/rs
	return (1-Flens(r,rs))*(8/(y*y*y) - 20/y + 15*y)
end

function gammaLens(r,rs)
	y = r/rs
	check = y-1
	check2 = y*y - 1
	if check < 0 && check2 < 0
		return 8*atanh(sqrt((1-y)/(1+y)))/(y*y*sqrt(1-y*y)) + 4*log(y/2)/(y*y) - 2/(y*y-1) + 4*atanh(sqrt((1-y)/(1+y)))/((y*y-1)*sqrt(1-y*y))
	elseif check > 0 && check2 > 0
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

function kappafunc_sis(theta, theta_e, dl)
	theta = theta / dl
	return  theta_e ./ (2 .* theta) 
end

function gammafunc_nfw(r,rs,gs)
	rs = rs*u"kpc"
	invGs = (1/gs)*u"kpc"
	return rs ./ invGs .* gammaLens.(r,rs)
end

function gammafunc_sis(theta, theta_e, dl)
	theta = theta / dl
	return  theta_e ./ (2 .* theta) 
end

function flexF_nfw(r,rs,fs)
	rs = rs*u"kpc"
	y = r/rs
	return (2 .* y .* Flens.(r,rs) .- Hlens.(r,rs)).* 2*fs ./ (y .^2 .-1).^2
end

function flexF_sis(theta, theta_e, dl)
	theta = theta / dl
	return  theta_e ./ (2 .* theta .^2)
end

function flexG_nfw(r,rs,fs) 
	rs = rs*u"kpc"
	y = r/rs
	return 2*fs .* (( log.(y/2) .* 8 ./ (y .* y .* y)) .+ ((3 ./ y) .* (1 .- 2 .* y .* y) .+ Glens.(r,rs)) ./ ((y .* y .- 1) .^2))
end

function flexG_sis(theta, theta_e, dl)
	theta = theta / dl
	return  3 .* theta_e ./ (2 .* theta .^2)
end


####################################################


####################################################



specbase = "alpha2jac_final"

arcm = 50f0
arcmdoub = 50

# 1 for jackknife, 0 for bootstrap
error_estimator = 1

# 2 for most massive galaxies, 1 for intermediate mass galaxies, 0 for low mass galaxies
massindic = 2

# 0 for central galaxies, 1 for satellite galaxies
satellite = 1

# NFW or SIS
modelHalo = "all"

ids = Int64[50,100,150,200,250,300,350,400,417,450]
zlist=Float64[0.21253133,0.3513419,0.5109155,0.7274435,1.0156995,1.4190572,1.8170974,2.4618596,2.79021206,3.913962228095617]
massivecentralmeanmass= Float64[11.638192255581975,12.962736628928043,13.141039914694,13.321059885996291,12.964841227475013,12.893973903248861,13.060299713205156,12.977746136120652,13.010365227427139,13.107033473239703]
leastmassivecentralmeanmass=Float64[10.785694241219115,10.43326496220352, 10.348051302604519,  10.279011566906112, 10.253380375911707, 10.242003300588081, 10.237021396526943,10.233305207105662, 10.232710124612193, 10.236236608629655] #10¹¹ Msun ?
massivesatellitemeanmass=Float64[10.093785133867357, 13.852519293406814, 13.848878736512019, 13.852062347990882, 13.512853174144187, 13.710894954066804, 13.922467501294063, 13.363746084518537, 13.416934113653896, 13.882345756005504]
leastmassivesatellitemeanmass=Float64[ 9.922009770623127, 3.0, 3.0, 3.0,  3.0, 3.0, 3.0,  3.0,  3.0, 3.0]

massivecentralmeanmass_subhalo = Float64[11.638192254085315, 12.962736628800336, 13.14103991463457, 13.32105988596674, 12.964841227425955, 12.893973903189963, 13.060299713165275, 12.977746136075575, 13.010365227385682, 13.107033473204401]
leastmassivecentralmeanmass_subhalo = Float64[10.785694233378985, 10.433264947421959, 10.348051284902992, 10.279011546303217, 10.253380354088051, 10.242003278193694, 10.237021373877202, 10.233305184262772, 10.23271010173815, 10.236236585939054]
massivesatellitemeanmass_subhalo = Float64[10.814911758074588, 10.717006737943423, 10.717910840074952, 10.723752635921732, 10.72153775407181, 10.812847733094795, 10.776229032950743, 10.7814438831561, 10.790639776654503, 10.777761971164157]
leastmassivesatellitemeanmass_subhalo = Float64[10.80813110659669, 10.541833642166472, 10.484367287226567, 10.458994770720563, 10.727680719783976, 10.730572389424118, 10.71421494527932, 10.84207383074002, 10.699474515356924, 10.766905483472398]

yscale = :log10
kwargs = Dict(:xlabel => "Distance", 
              :xscale => :log10, 
              :yscale => yscale, 
              :legend => yscale == :log10 ? :bottomleft : :topright,
			  :linewidth => 2,
			  :titlefontsize => 10,)

			  
kappa_fit_params_nfw = OrderedDict{Int64, Tuple{Vector{Float64}, Vector{Float64}}}()
gamma_fit_params_nfw = OrderedDict{Int64, Tuple{Vector{Float64}, Vector{Float64}}}()
F_fit_params_nfw = OrderedDict{Int64, Tuple{Vector{Float64}, Vector{Float64}}}()
G_fit_params_nfw = OrderedDict{Int64, Tuple{Vector{Float64}, Vector{Float64}}}()
m_c_dict_kappa = OrderedDict{Int64, Vector{Float64}}()
m_c_dict_gamma = OrderedDict{Int64, Vector{Float64}}()
m_c_dict_F = OrderedDict{Int64, Vector{Float64}}()
m_c_dict_G = OrderedDict{Int64, Vector{Float64}}()

kappa_fit_params_sis = OrderedDict{Int64, Tuple{Vector{Float64}, Vector{Float64}}}()
gamma_fit_params_sis = OrderedDict{Int64, Tuple{Vector{Float64}, Vector{Float64}}}()
F_fit_params_sis = OrderedDict{Int64, Tuple{Vector{Float64}, Vector{Float64}}}()
G_fit_params_sis = OrderedDict{Int64, Tuple{Vector{Float64}, Vector{Float64}}}()

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


if satellite == 0 && massindic == 2
	meanmass = massivecentralmeanmass_subhalo
elseif satellite == 0 && massindic == 0
	meanmass = leastmassivecentralmeanmass_subhalo
elseif satellite == 1 && massindic == 2
	meanmass = massivesatellitemeanmass_subhalo
elseif satellite == 1 && massindic == 0
	meanmass = leastmassivesatellitemeanmass_subhalo
end	


@info @sprintf("Starting analysis for the %s %s galaxies",masstitle,galtitle)
for (i, val) in enumerate(ids)
	model = modelHalo
    z_now = zlist[i]
	meanmass_id = meanmass[i]*cosmo.h*1e11
	
	specifier = specbase

	half_d = 0.5*angular_diameter_dist(cosmo, z_now)
	
	xi(z) = angular_diameter_dist(cosmo,z) - half_d
	
	z_half = find_zero(xi,z_now/i)
	comov_d = angular_diameter_dist(cosmo, z_half)
	plot_z_now = @sprintf("%.3f", z_now)
	plot_z_half = @sprintf("%.3f", z_half)

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

	if model == "NFW"
		specifier = "NFW_"*specifier
	elseif model == "SIS"
		specifier = "SIS_"*specifier
	else model == "all"
		specifier = "NFW-SIS_"*specifier
	end

	specifier = "noSimMass_"*specifier

	G_temp = Sres[1,:,:] ./ Sres[5,:,:]
	F_temp = Sres[2,:,:] ./ Sres[5,:,:]
	kappa_temp = Sres[3,:,:] ./ Sres[5,:,:]
	gamma_temp = Sres[4,:,:] ./ Sres[5,:,:]
	
	G = mean(G_temp, dims=2) 
	F = mean(F_temp, dims=2)
	gamma = mean(gamma_temp, dims=2)
	kappa = mean(kappa_temp, dims=2)
	
	if G[1] < 0 
		G = -G										
	end

	if F[1] < 0 
		F = -F
	end

	if gamma[1] < 0 
		gamma = -gamma
	end

	if kappa[1] < 0 
		kappa = -kappa
	end

	"""nan_indices_F = findall(vec(isnan.(F) .& F .< 0))	
	nan_indices_G = findall(vec(isnan.(G) .& G .< 0))	
	nan_indices_gamma = findall(vec(isnan.(gamma) .& gamma .< 0))	
	nan_indices_kappa = findall(vec(isnan.(kappa) .& kappa .< 0))	"""
	rbin_centers_F = rbin_centers[vec(isfinite.(F) .& (F .> 0))]
	rbin_centers_G = rbin_centers[vec(isfinite.(G) .& (G .> 0))]
	rbin_centers_gamma = rbin_centers[vec(isfinite.(gamma) .& (gamma .> 0))]
	rbin_centers_kappa = rbin_centers[vec(isfinite.(kappa) .& (kappa .> 0))]

	F = F[vec(isfinite.(F) .& (F .> 0))]
	G = G[vec(isfinite.(G) .& (G .> 0))]
	gamma = gamma[vec(isfinite.(gamma) .& (gamma .> 0))]
	kappa = kappa[vec(isfinite.(kappa) .& (kappa .> 0))]


	if error_estimator == 1 							#JACKKNIFE 
		#specifier = "jackknifeerrors_"*specifier

		n_jackknife = 5000
		n_exclude = floor(size(Sres, 3) ÷ 100)
		G_jackknife = zeros(size(G_temp, 1), n_jackknife)
		F_jackknife = zeros(size(F_temp, 1), n_jackknife)
		gamma_jackknife = zeros(size(gamma_temp, 1), n_jackknife)
		kappa_jackknife = zeros(size(kappa_temp, 1), n_jackknife)


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

		for i in 1:n_jackknife
		    for j in 1:size(G_jackknife, 1)
			G_jackknife_i[j, i] = (G_jackknife[j, i] - mean_Gjk[j])^2
		    end
		    for j in 1:size(F_jackknife, 1)
			F_jackknife_i[j, i] = (F_jackknife[j, i] - mean_Fjk[j])^2
		    end
		    for j in 1:size(gamma_jackknife, 1)
			gamma_jackknife_i[j, i] = (gamma_jackknife[j, i] - mean_gammajk[j])^2
                    end
		    for j in 1:size(kappa_jackknife, 1)
			kappa_jackknife_i[j, i] = (kappa_jackknife[j, i] - mean_kappajk[j])^2
	 	    end
		end

		G_error = sqrt.((n_jackknife - 1) * mean(G_jackknife_i, dims=2))
		F_error = sqrt.((n_jackknife - 1) * mean(F_jackknife_i, dims=2))
		gamma_error = sqrt.((n_jackknife - 1) * mean(gamma_jackknife_i, dims=2))
		kappa_error = sqrt.((n_jackknife - 1) * mean(kappa_jackknife_i, dims=2))
	
	end

	#if get(kwargs, :yscale, :identity) == :log10
	#	@info("Absolute value for log")
	#	G = abs.(G)
	#	F = abs.(F)
	#	gamma = abs.(gamma)
	#	kappa = abs.(kappa)
	#end
	
	#G = abs.(G)
	#F = abs.(F)
	#gamma = abs.(gamma)
	#kappa = abs.(kappa)

	rbin_centers_G = comov_d*rbin_centers_G*3.14159/(180)
	rbin_centers_G = uconvert.(u"kpc", rbin_centers_G) 

	rbin_centers_F = comov_d*rbin_centers_F*3.14159/(180) 
	rbin_centers_F = uconvert.(u"kpc", rbin_centers_F) 

	rbin_centers_gamma = comov_d*rbin_centers_gamma*3.14159/(180) 
	rbin_centers_gamma = uconvert.(u"kpc", rbin_centers_gamma) 

	rbin_centers_kappa = comov_d*rbin_centers_kappa*3.14159/(180) 
	rbin_centers_kappa = uconvert.(u"kpc", rbin_centers_kappa) 


	pG = scatter(rbin_centers_G, G, yerr = G_error, ylabel="G",  title = "\nG for ID = $val "*L"$z_{plane} = $"*"$plot_z_now "*L"$z_{deflectors} = $"*"$plot_z_half"; kwargs...)
	pF = scatter(rbin_centers_F, F, yerr = F_error, ylabel="F",  title = "\nF fit for ID = $val "*L"$z_{plane} = $"*"$plot_z_now "*L"$z_{deflectors} = $"*"$plot_z_half"; kwargs...)
	pgamma = scatter(rbin_centers_gamma, gamma, yerr = gamma_error, ylabel=L"\gamma",  title = "\nShear for ID = $val "*L"$z_{plane} = $"*"$plot_z_now "*L"$z_{deflectors} = $"*"$plot_z_half"; kwargs...)
	pkappa = scatter(rbin_centers_kappa, kappa, yerr = kappa_error, ylabel=L"\kappa",  title = "\nConvergence fit for ID = $val "*L"$z_{plane} = $"*"$plot_z_now "*L"$z_{deflectors} = $"*"$plot_z_half"; kwargs...)

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
	
	pall = scatter(rbin_centers_G, G, yerr = G_error, ylabel="Value", title="correlations for $val", label=L"G_+"; kwargs... )
	scatter!(rbin_centers_F, F, yerr = F_error, label=L"F_+"; kwargs... )
	scatter!(rbin_centers_gamma, gamma, yerr = gamma_error, label=L"\gamma"; kwargs... )
	scatter!(rbin_centers_kappa, kappa, yerr = kappa_error, label=L"\kappa"; kwargs... )

	savefig(pall,allpath)


	
	modelkappa_nfw = (x, p) -> kappafunc_nfw(x, p[1], p[2])
	modelgamma_nfw = (x, p) -> gammafunc_nfw(x, p[1], p[2])
	modelF_nfw = (x, p) -> flexF_nfw(x, p[1], p[2])
	modelG_nfw = (x, p) -> flexG_nfw(x, p[1], p[2])
	p0_nfw = [1, 1]

	modelkappa_sis = (x, p) -> kappafunc_sis(x, p[1], comov_d)
	modelgamma_sis = (x, p) -> gammafunc_sis(x, p[1], comov_d)
	modelF_sis = (x, p) -> flexF_sis(x, p[1], comov_d)
	modelG_sis = (x, p) -> flexG_sis(x, p[1], comov_d)
	p0_sis = [10]

	# For NFW models, if parameters need to be positive, transform them to log-space
	modelkappa_nfw = (x, log_p) -> kappafunc_nfw(x, exp(log_p[1]), exp(log_p[2]))
	modelgamma_nfw = (x, log_p) -> gammafunc_nfw(x, exp(log_p[1]), exp(log_p[2]))
	modelF_nfw = (x, log_p) -> flexF_nfw(x, exp(log_p[1]), exp(log_p[2]))
	modelG_nfw = (x, log_p) -> flexG_nfw(x, exp(log_p[1]), exp(log_p[2]))
	p0_nfw = log.(p0_nfw)
	
	# For SIS models, if parameters need to be positive, transform them similarly
	modelkappa_sis = (x, log_p) -> kappafunc_sis(x, exp(log_p[1]), comov_d)
	modelgamma_sis = (x, log_p) -> gammafunc_sis(x, exp(log_p[1]), comov_d)
	modelF_sis = (x, log_p) -> flexF_sis(x, exp(log_p[1]), comov_d)
	modelG_sis = (x, log_p) -> flexG_sis(x, exp(log_p[1]), comov_d)
	p0_sis = log.(p0_sis)


	fit_cutoff = 5e3*u"kpc"
	fit_rbin_centers_kappa = rbin_centers_kappa[(rbin_centers_kappa .< fit_cutoff)]
	fit_rbin_centers_gamma = rbin_centers_gamma[(rbin_centers_gamma .< fit_cutoff)]
	fit_rbin_centers_F = rbin_centers_F[(rbin_centers_F .< fit_cutoff)]
	fit_rbin_centers_G = rbin_centers_G[(rbin_centers_G .< fit_cutoff)]

	fit_kappa = kappa[(rbin_centers_kappa .< fit_cutoff)]
	fit_gamma = gamma[(rbin_centers_gamma .< fit_cutoff)]
	fit_F = F[(rbin_centers_F .< fit_cutoff)]
	fit_G = G[(rbin_centers_G .< fit_cutoff)]

	#Fitting lensing quantities
	kappafitpath = @sprintf("corr_plots/fits/kappa/kappa_fit_%s_%d_arcm_id_%.4d.png",specifier,arcmdoub,val)
	gammafitpath = @sprintf("corr_plots/fits/gamma/gamma_fit_%s_%d_arcm_id_%.4d.png",specifier,arcmdoub,val)
	Ffitpath = @sprintf("corr_plots/fits/F/F_fit_%s_%d_arcm_id_%.4d.png",specifier,arcmdoub,val)
	Gfitpath = @sprintf("corr_plots/fits/G/G_fit_%s_%d_arcm_id_%.4d.png",specifier,arcmdoub,val)

	kappaFit_nfw = curve_fit(modelkappa_nfw, fit_rbin_centers_kappa, fit_kappa, p0_nfw)
	#kappaFitErrors_nfw = estimate_errors(kappaFit_nfw, 0.95)
	kappaFitErrors_nfw = [0,0]

	gammaFit_nfw = curve_fit(modelgamma_nfw, fit_rbin_centers_gamma, fit_gamma, p0_nfw)
	#gammaFitErrors_nfw = estimate_errors(gammaFit_nfw, 0.95)
	gammaFitErrors_nfw = [0,0]

	FFit_nfw = curve_fit(modelF_nfw, fit_rbin_centers_F, fit_F, p0_nfw)
	#FFitErrors_nfw = estimate_errors(FFit_nfw, 0.95)
	FFitErrors_nfw = [0,0]

	GFit_nfw = curve_fit(modelG_nfw, fit_rbin_centers_G, fit_G, p0_nfw)
	#GFitErrors_nfw = estimate_errors(GFit_nfw, 0.95)
	GFitErrors_nfw = [0,0]

	kappaFit_sis = curve_fit(modelkappa_sis, fit_rbin_centers_kappa, fit_kappa, p0_sis)
	kappaFitErrors_sis = estimate_errors(kappaFit_sis, 0.95)

	gammaFit_sis = curve_fit(modelgamma_sis, fit_rbin_centers_gamma, fit_gamma, p0_sis)
	gammaFitErrors_sis = estimate_errors(gammaFit_sis, 0.95)

	FFit_sis = curve_fit(modelF_sis, fit_rbin_centers_F, fit_F, p0_sis)
	FFitErrors_sis = estimate_errors(FFit_sis, 0.95)

	GFit_sis = curve_fit(modelG_sis, fit_rbin_centers_G, fit_G, p0_sis)
	GFitErrors_sis = estimate_errors(GFit_sis, 0.95)

	rvec_kappa = range(ustrip(rbin_centers_kappa[1]), stop=ustrip(rbin_centers_kappa[end]) ,length=1000)*u"kpc" 
	rvec_gamma = range(ustrip(rbin_centers_gamma[1]), stop=ustrip(rbin_centers_gamma[end]) ,length=1000)*u"kpc"
	rvec_F = range(ustrip(rbin_centers_F[1]), stop=ustrip(rbin_centers_F[end]) ,length=1000)*u"kpc"
	rvec_G = range(ustrip(rbin_centers_G[1]), stop=ustrip(rbin_centers_G[end]) ,length=1000)*u"kpc"

	"""fittedKappa_nfw = modelkappa_nfw(rvec_kappa, collect(kappaFit_nfw.param))
	fittedGamma_nfw = modelgamma_nfw(rvec_gamma, collect(gammaFit_nfw.param))
	fittedF_nfw = modelF_nfw(rvec_F, collect(FFit_nfw.param))
	fittedG_nfw = modelG_nfw(rvec_G, collect(GFit_nfw.param))

	fittedKappa_sis = modelkappa_sis(rvec_kappa, collect(kappaFit_sis.param))
	fittedGamma_sis = modelgamma_sis(rvec_gamma, collect(gammaFit_sis.param))
	fittedF_sis = modelF_sis(rvec_F, collect(FFit_sis.param))
	fittedG_sis = modelG_sis(rvec_G, collect(GFit_sis.param))"""

	fittedKappa_nfw = modelkappa_nfw(rvec_kappa, kappaFit_nfw.param)
	fittedGamma_nfw = modelgamma_nfw(rvec_gamma, gammaFit_nfw.param)
	fittedF_nfw = modelF_nfw(rvec_F, FFit_nfw.param)
	fittedG_nfw = modelG_nfw(rvec_G, GFit_nfw.param)

	fittedKappa_sis = modelkappa_sis(rvec_kappa, kappaFit_sis.param)
	fittedGamma_sis = modelgamma_sis(rvec_gamma, gammaFit_sis.param)
	fittedF_sis = modelF_sis(rvec_F, FFit_sis.param)
	fittedG_sis = modelG_sis(rvec_G, GFit_sis.param)


	(mass, c) = ksrs2mc(exp.(kappaFit_nfw.param[1]), exp.(kappaFit_nfw.param[2]), 1.0, runit = :kpc, z=z_half)
	m_c_dict_kappa[val] = [mass, c]
	
	(mass, c) = ksrs2mc(exp.(gammaFit_nfw.param[1]), exp.(gammaFit_nfw.param[1]), 1.0, runit = :kpc, z=z_half)
	m_c_dict_gamma[val] = [mass, c]

	(mass, c) = ksrs2mc(exp.(FFit_nfw.param[1]), exp.(FFit_nfw.param[2]), 1.0, runit = :kpc, z=z_half)
	m_c_dict_F[val] = [mass, c]

	(mass, c) = ksrs2mc(exp.(GFit_nfw.param[1]), exp.(GFit_nfw.param[2]), 1.0, runit = :kpc, z=z_half)
	m_c_dict_G[val] = [mass, c]

	
	(ks, rs) = m2ksrs(cm_itp, meanmass_id, z_half)
	rs = ustrip(ustrip(rs)*uconvert(u"kpc",comov_d)*3.14159/(180))
	KapGamSimparam = log.([rs, ustrip(ks)/2])
	print("meanmass_id = ", meanmass_id, "\n")
	print("rs = ", rs, "\n")
	print("ks = ", ks, "\n")

	fs = ks*ustrip(comov_d)/rs
	FGSimparam = log.([rs, fs])
	print("fs = ", fs, "\n")
	

	"""kappa_fit_params_nfw[val] = (kappaFit_nfw.param, kappaFitErrors_nfw)
	gamma_fit_params_nfw[val] = (gammaFit_nfw.param, gammaFitErrors_nfw)
	F_fit_params_nfw[val] = (FFit_nfw.param, FFitErrors_nfw)
	G_fit_params_nfw[val] = (GFit_nfw.param, GFitErrors_nfw)

	kappa_fit_params_sis[val] = (kappaFit_sis.param, kappaFitErrors_sis)
	gamma_fit_params_sis[val] = (gammaFit_sis.param, gammaFitErrors_sis)
	F_fit_params_sis[val] = (FFit_sis.param, FFitErrors_sis)
	G_fit_params_sis[val] = (GFit_sis.param, GFitErrors_sis)"""

	kappa_fit_params_nfw[val] = (exp.(kappaFit_nfw.param), exp.(kappaFitErrors_nfw))
	gamma_fit_params_nfw[val] = (exp.(gammaFit_nfw.param), exp.(gammaFitErrors_nfw))
	F_fit_params_nfw[val] = (exp.(FFit_nfw.param), exp.(FFitErrors_nfw))
	G_fit_params_nfw[val] = (exp.(GFit_nfw.param), exp.(GFitErrors_nfw))

	kappa_fit_params_sis[val] = (exp.(kappaFit_sis.param), exp.(kappaFitErrors_sis))
	gamma_fit_params_sis[val] = (exp.(gammaFit_sis.param), exp.(gammaFitErrors_sis))
	F_fit_params_sis[val] = (exp.(FFit_sis.param), exp.(FFitErrors_sis))
	G_fit_params_sis[val] = (exp.(GFit_sis.param), exp.(GFitErrors_sis))
	#if get(kwargs, :yscale, :identity) == :log10
	#	fittedKappa = abs.(fittedKappa)
	#end
	
	

	k = scatter(rbin_centers_kappa, kappa,yerr = kappa_error, label = L"\bar\kappa", ylabel = L"\kappa", title = "\nConvergence fit for ID = $val "*L"$z_{plane} = $"*"$plot_z_now "*L"$z_{deflectors} = $"*"$plot_z_half"; kwargs...)
	#plot!(rvec_kappa, modelkappa_nfw(rvec_kappa,KapGamSimparam), label="Simulation mass NFW "*L"\kappa"*" profile"; kwargs...)
	plot!(rvec_kappa, fittedKappa_nfw, label="NFW "*L"\kappa"*" profile fit"; kwargs...)
	plot!(rvec_kappa, fittedKappa_sis, label="SIS "*L"\kappa"*" profile fit"; kwargs...)
	savefig(k,kappafitpath)

	gam = scatter(rbin_centers_gamma, gamma,yerr = gamma_error, label = L"\bar\gamma", ylabel = L"\gamma",  title = "\nShear fit for ID = $val "*L"$z_{plane} = $"*"$plot_z_now "*L"$z_{deflectors} = $"*"$plot_z_half"; kwargs...)
	#plot!(rvec_gamma, modelgamma_nfw(rvec_gamma,KapGamSimparam), label="Simulation mass NFW "*L"\gamma"*" profile"; kwargs...)
	plot!(rvec_gamma, fittedGamma_nfw, label="NFW "*L"\gamma"*" profile fit"; kwargs...)
	plot!(rvec_gamma, fittedGamma_sis, label="SIS "*L"\gamma"*" profile fit"; kwargs...)
	savefig(gam,gammafitpath)

	f = scatter(rbin_centers_F, F,yerr = F_error, label = L"\bar{F}", ylabel = "F in 1/arcsecond", title = "\nF fit for ID = $val "*L"$z_{plane} = $"*"$plot_z_now "*L"$z_{deflectors} = $"*"$plot_z_half"; kwargs...)
	#plot!(rvec_F, modelF_nfw(rvec_F,FGSimparam), label="Simulation mass NFW F profile"; kwargs...)
	plot!(rvec_F, fittedF_nfw, label="NFW F profile fit"; kwargs...)
	plot!(rvec_F, fittedF_sis, label="SIS F profile fit"; kwargs...)
	savefig(f,Ffitpath)

	g = scatter(rbin_centers_G, G,yerr = G_error, label = L"\bar{G}", ylabel = "G in 1/kpc",title = "\nG fit for ID = $val "*L"$z_{plane} = $"*"$plot_z_now "*L"$z_{deflectors} = $"*"$plot_z_half"; kwargs...)
	#plot!(rvec_G, modelG_nfw(rvec_G,FGSimparam), label="Simulation mass NFW G profile"; kwargs...)
	plot!(rvec_G, fittedG_nfw, label="NFW G profile fit"; kwargs...)
	plot!(rvec_G, fittedG_sis, label="SIS G profile fit"; kwargs...)
	savefig(g,Gfitpath)
end


@info "Kappa SIS parameters"
println(kappa_fit_params_sis)
@info "Gamma SIS parameters"
println(gamma_fit_params_sis)
@info "F SIS parameters"
println(F_fit_params_sis)
@info "G SIS parameters"
println(G_fit_params_sis)

@info "Kappa NFW parameters"
println(kappa_fit_params_nfw)
@info "Gamma NFW parameters"
println(gamma_fit_params_nfw)
@info "F NFW parameters"
println(F_fit_params_nfw)
@info "G NFW parameters"
println(G_fit_params_nfw)

@info "Mass and concentration values for kappa"
println(m_c_dict_kappa)
@info "Mass and concentration values for gamma"
println(m_c_dict_gamma)
@info "Mass and concentration values for F"
println(m_c_dict_F)
@info "Mass and concentration values for G"
println(m_c_dict_G)


cvec = range(0.1, stop=2.0, length=100)
δ_c_vec = δ_c.(cvec)
masses = []
c_values = []
for (key, value) in m_c_dict_kappa
    push!(masses, value[1])
    push!(c_values, value[2])
end
#plot(cvec, δ_c_vec, xlabel="c", ylabel=L"\delta_c", title=L"\delta_c"*" vs c", label=L"\delta_c vs c", legend=true; kwargs...)
scatter(c_values, masses, label="Mass vs c", legend=true; kwargs...)


