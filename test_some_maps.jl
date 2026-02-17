include("raytrace_tk.jl")

#outdir="/scratchb08/gavazzi/lensing_maps"
#ids = reverse(20:20:500)
ids = Int64[50,100,150,200,250,300,350,400,417,450]
totest= String["rot","a1crossderivative","a2crossderivative"]

for id in ids
	for test in totest
	
	    outdir="/net/GECO/nas12c/euclid/TestWorkFolder/test_maps/"*test
	    isfile(outdir) || mkpath(outdir)
	    dx=1.0f0
	    mode="OBB"
	    rf=5
	    update=false

       	    tbl_io = open(joinpath(outdir,"id_z_mapping.csv"),"w")
       	    println(tbl_io,"id,redshift")

	
	    fn = full_bin_name( id, mode=mode, size=dx )

	    if ( !isfile(fn) & !update )
		@warn "$fn doesn't exist, skipping."
		continue
	    end

	    startoutfilenametemp = "HAGN-"*test*"-test"
	    endoutfilenametemp = @sprintf("_%04d.fits",id)
	    ofnp = startoutfilenametemp*endoutfilenametemp
	    ofn = joinpath( outdir, ofnp )

	    if ( isfile(ofn) & !update )
		@warn "$ofn already exists. Skipping!"
		continue
	    end

	    ######### Deflection and derived information
	    @time (alpha,zz) = read_bin( fn )
	    @printf( tbl_io, "%d,%0.3f\n", id, zz )
	    nn=size(alpha)
	    @info "Done reading deflection map for plane $id"
	    @time jac=alpha2jac( alpha, scale=dx )
	    @info "Done computing first order derivatives !"
	    @time jac=jac2der( jac )
	    @info "Done computing full jacobian for flexion !"
	    @time smu = alpha_mu_2_mus( alpha, dx ./ collect(size(alpha)) )
	    alpha=nothing ; GC.gc(true)	
	    @info "Done computing jacobian matrix map"
	    if rf>1
		@time begin
		    #kappa = rebin( jac2kappa(jac), factor=rf )
		    #gamma1 = rebin( jac2gamma1(jac), factor=rf )
		    #gamma2 = rebin( jac2gamma2(jac), factor=rf )
		    #f1flexion = rebin(jac2F1(jac), factor=rf)
		    #gflexion = rebin(jac2G(jac), factor=rf)
		if test == "rot"
		    rot = rebin( jac2rot(jac), factor=rf )
		elseif test == "a1crossderivative"
		    a1CD = rebin( jac2crossderivativea1(jac), factor=rf )
		elseif test == "a2crossderivative"
		    a2CD = rebin( jac2crossderivativea2(jac), factor=rf )
		end
		    #imu = rebin( jac2muinv(jac), factor=rf )
		    #smu = rebin( smu, factor=rf )
		end
	    else
		@time begin
		    #kappa = jac2kappa(jac)
		    #gamma1 = jac2gamma1(jac)
		    #gamma2 = jac2gamma2(jac)
		    #rot = jac2rot(jac)
		    #imu = jac2muinv(jac)
		if test == "rot"
		    rot = jac2rot(jac)
		elseif test == "a1crossderivative"
		    a1CD = jac2crossderivativea1(jac)
		elseif test == "a2crossderivative"
		    a2CD = jac2crossderivativea2(jac)
		end
		    #f1flexion = jac2F1(jac)
		    #gflexion = jac2G(jac)
		    ## smu is ok
		end
	    end
	    jac=nothing ; GC.gc(true)
	    @info "Done computing rebinned flexion maps"

	    nn = ( nn[1]/rf, nn[2]/rf )
	    @info "Storing in file $ofn"

	    FITS( ofn, "w" ) do io
		hd = wcs_header( zz, dx, nn )
		#hd["MAP"] = "kappa"
		#write( io, kappa, header=hd )
		#hd["MAP"] = "gamma1"
		#write( io, gamma1, header=hd )
		#hd["MAP"] = "gamma2"
		#write( io, gamma2, header=hd )
		#hd["MAP"] = "F1"
		#write( io, f1flexion, header=hd )
		if test == "rot"
			hd["MAP"] = "rot"
			write( io, rot, header=hd )
		elseif test == "a1crossderivative"
			hd["MAP"] = "A1CROSSDER"
			write( io, a1CD, header=hd )
		elseif test == "a2crossderivative"
			hd["MAP"] = "A2CROSSDER"
			write( io, a2CD, header=hd )
		end				
		#hd["MAP"] = "G"
		#write( io, gflexion, header=hd )
		#hd["MAP"] = "mu_img"
		#write( io, imu, header=hd )
		#hd["MAP"] = "mu_src"
		#write( io, smu, header=hd )
	    end
		close(tbl_io)
	end


end
