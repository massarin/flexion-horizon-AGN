using JLD2, LaTeXStrings, Cosmology, Plots, LsqFit
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


ids = Int64[50,100,150,200,250,300,350,400,417,450]
zlist=Float64[0.21253133,0.3513419,0.5109155,0.7274435,1.0156995,1.4190572,1.8170974,2.4618596,2.79021206,3.913962228095617]


#ids=reverse(ids)
#zlist=reverse(zlist)

arcm = 5f0
arcmdoub = 5

####################################################
function atanhf(x)
	return 1-2*atanh(sqrt((1-x)/(1+x)))/sqrt(1-x*x)
end


function atanf(x)
	return 1-2*atan(sqrt((x-1)/(1+x)))/sqrt(x*x-1)
end


function ffunc(x)
	if x > 1.0
		return atanf(x)
	elseif x < 1.0
		return atanhf(x)
	else 
		return 0
	end
end


function hfunc(x)
	return 1-x*ffunc(x)-1/x
end 


function gfunc(x)
	return (1-ffunc(x))*(8/(x*x*x) - 20/x + 15*x)
end
####################################################

function kappafunc(x, rs, ks)
	return 2 .* ks .* ffunc.(x ./ rs) ./ ((x .* x) ./rs .- 1)
end

function flexF(x, Fs)
	return -(2 .* x .* ffunc.(x) .- hfunc.(x)) .* 2 .* Fs ./ (x .* x .- 1).^ 2
end

function flexG(x, Fs) 
	return 2*Fs .* (log.(x ./ 2) .* 8 ./ (x .* x .* x) .+ ((3 ./ x) .* (1 .- 2 .* x .* x) .+ gfunc.(x)) ./ (x .* x .- 1) .^ 2)
end

####################################################

####################################################

model(x,p) = kappafunc(x,p[1],p[2])
p0 = [1.0,1.0]

ids = Int64[50,100,150,200,250,300,350,400,417,450]
ids = Int64[50]

fit_params = Dict{Int64, Vector{Float64}}()

p = scatter([],[],xscale=:log10,yscale=:log10,legend=false)
for val in ids

	datafile = @sprintf("Data/corr_data/Sres_rbin_alpha2jacwcrossterms_id_0%d_arcm%d.jld2",val,arcmdoub)
	@load datafile Sres rbin 
	@eval $(Symbol("Sres$val")) = $Sres
	@eval $(Symbol("rbin$val")) = $rbin
	rbin_centers = 0.5*(rbin[1:end-1]+rbin[2:end])

	G_temp = Sres[1,:,:] ./ Sres[5,:,:]
	F_temp = Sres[2,:,:] ./ Sres[5,:,:]
	kappa_temp = Sres[3,:,:] ./ Sres[5,:,:]
	gamma_temp = Sres[4,:,:] ./ Sres[5,:,:]

	G = mean(G_temp, dims=2)
	F = mean(F_temp, dims=2)
	gamma = mean(gamma_temp, dims=2)
	kappas = vec(mean(kappa_temp, dims=2))
	
	
	kappaFit = curve_fit(model, rbin_centers, kappas, p0)
	
	fit_params[val]=kappaFit.param

	fittedKappa = model(rbin_centers,kappaFit.param)
	scatter!(p,rbin_centers, abs.(kappas),xscale=:log10,yscale=:log10)
	plot!(p,rbin_centers, abs.(fittedKappa),xscale=:log10,yscale=:log10)
	
end
display(p)


