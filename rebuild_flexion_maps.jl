include("raytrace_tk.jl")

#outdir="/scratchb08/gavazzi/lensing_maps"
outdir="/net/GECO/nas12c/euclid/TestWorkFolder/lensing_maps/flexion_maps"
isfile(outdir) || mkpath(outdir)
dx=1.0f0
mode="OBB"
rf=5
update=false

tbl_io = open(joinpath(outdir,"id_z_mapping.csv"),"w")
println(tbl_io,"id,redshift")
ids = reverse(20:20:500)
ids = Int64[50,100,150,200,250,300,350,400,417,450]
for id in ids
    fn = full_bin_name( id, mode=mode, size=dx )
    if ( !isfile(fn) & !update )
	@warn "$fn doesn't exist, skipping."
	continue
    end
    ofn = joinpath( outdir,
                    @sprintf("HAGN-flexion-lightcone_%04d.fits",id) )
    if ( isfile(ofn) & !update )
        @warn "$ofn already exists. Skipping!"
        continue
    end
    ######### Deflection and derived information
    @time (alpha,zz) = read_bin( fn )
    @printf( tbl_io, "%d,%0.3f\n", id, zz )
    nn=size(alpha)
    @info "Done reading deflection map for plane $id"
    @time jac=alpha2jac_flexion( alpha, scale=dx )
    @info "Computing source plane total magnification"
    @time smu = alpha_mu_2_mus( alpha, dx ./ collect(size(alpha)) )
    alpha=nothing ; GC.gc(true)
    @info "Done computing jacobian matrix map"
    if rf>1
        @time begin
            #kappa = rebin( jac2kappa(jac), factor=rf )
            #gamma1 = rebin( jac2gamma1(jac), factor=rf )
            #gamma2 = rebin( jac2gamma2(jac), factor=rf )
	    fflexion = rebin(jac2F(jac), factor=rf)
	    gflexion = rebin(jac2G(jac), factor=rf)
            #rot = rebin( jac2rot(jac), factor=rf )
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
	    fflexion = jac2F(jac)
	    gflexion = jac2G(jac)
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
	hd["MAP"] = "F"
        write( io, fflexion, header=hd )
	hd["MAP"] = "G"
        write( io, gflexion, header=hd )
        #hd["MAP"] = "rot"
        #write( io, rot, header=hd )
        #hd["MAP"] = "mu_img"
	#write( io, imu, header=hd )
        #hd["MAP"] = "mu_src"
        #write( io, smu, header=hd )
    end

end

close(tbl_io)
