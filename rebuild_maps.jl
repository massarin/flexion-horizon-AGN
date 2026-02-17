include("raytrace_tk.jl")

#outdir="/scratchb08/gavazzi/lensing_maps"
outdir="/net/GECO/nas12c/euclid/TestWorkFolder/lensing_maps/all"
isfile(outdir) || mkpath(outdir)
dx=1.0f0
mode="OBB"
rf=5
update=true

tbl_io = open(joinpath(outdir,"id_z_mapping.csv"),"w")
println(tbl_io,"id,redshift")
ids = reverse(20:20:500)
ids = Int64[50,100,150,200,250,300,350,400,417,450]

Threads.@threads for id in ids
    fn = full_bin_name( id, mode=mode, size=dx )
    if ( !isfile(fn) & !update )
	@warn "$fn doesn't exist, skipping."
	continue
    end
    ofn = joinpath( outdir,
                    @sprintf("HAGN-lensingderiv2ord_%04d.fits",id) )
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
    @info "Computing source plane total magnification"
    @time smu = alpha_mu_2_mus( alpha, dx ./ collect(size(alpha)) )
    alpha=nothing ; GC.gc(true)
    @info "Done computing jacobian matrix map"
    if rf>1
        @time begin
            kappa = rebin( jac2kappa(jac), factor=rf )
            gamma1 = rebin( jac2gamma1(jac), factor=rf )
            gamma2 = rebin( jac2gamma2(jac), factor=rf )
            F1 = rebin( jac2F1(jac), factor=rf )
            F2 = rebin( jac2F2(jac), factor=rf )
            G1 = rebin( jac2G1(jac), factor=rf )
            G2 = rebin( jac2G2(jac), factor=rf )
            rot = rebin( jac2rot(jac), factor=rf )
            a1crossder =  rebin(jac2crossderivativea1(jac), factor=rf)
            a2crossder =  rebin(jac2crossderivativea2(jac), factor=rf)
            imu = rebin( jac2muinv(jac), factor=rf )
            smu = rebin( smu, factor=rf )
        end
    else
        @time begin
            kappa = jac2kappa(jac)
            gamma1 = jac2gamma1(jac)
            gamma2 = jac2gamma2(jac)
	    F1 = jac2F1(jac)
            F2 = jac2F2(jac)
            G1 = jac2G1(jac)
            G2 = jac2G2(jac)
            rot = jac2rot(jac)
            a1crossder = jac2crossderivativea1(jac)
            a2crossder = jac2crossderivativea2(jac)
            imu = jac2muinv(jac)
            ## smu is ok
        end
    end
    jac=nothing ; GC.gc(true)
    @info "Done computing rebinned convergence & shear maps"

    nn = ( nn[1]/rf, nn[2]/rf )
    @info "Storing in file $ofn"

    FITS( ofn, "w" ) do io
        hd = wcs_header( zz, dx, nn )
        hd["KAPPA"] = "kappa"
        write( io, kappa, header=hd )
        hd["GAMMA_1"] = "gamma1"
        write( io, gamma1, header=hd )
        hd["GAMMA_2"] = "gamma2"
        write( io, gamma2, header=hd )
        hd["F1"] = "F1"
        write( io, F1, header=hd )
        hd["F2"] = "F2"
        write( io, F2, header=hd )
        hd["G1"] = "G1"
        write( io, G1, header=hd )
        hd["G2"] = "G2"
        write( io, G2, header=hd )
        hd["ROT"] = "rot"
        write( io, rot, header=hd )
        hd["A1CROSSDER"] = "a1crossder"
        write( io, a1crossder, header=hd )
        hd["A2CROSSDER"] = "a2crossder"
        write( io, a2crossder, header=hd )
        hd["mu_img"] = "mu_img"
	write( io, imu, header=hd )
        hd["mu_src"] = "mu_src"
        write( io, smu, header=hd )
    end
	kappa = gamma1 = gamma2 = F1 = F2 = G1 = G2 =nothing ; GC.gc(true)
end

close(tbl_io)
