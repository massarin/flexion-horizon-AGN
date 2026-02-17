using ProgressMeter, Distributed, ClusterManagers, FITSIO, Random, Dates

#(nprocs()<2) && addprocs(LocalAffinityManager(np=2,mode=BALANCED,affinities=[0,1]))

#nprocs()>1 && map( ip-> remotecall_fetch( ()->run(`taskset -p $(2^(ip-1)) $(getpid())`),workers()[ip]), 1:nworkers() )


@everywhere using Dates, Printf, FortranFiles, Interpolations


@everywhere begin

    const dx=1.0f0
    const IMADIR="/scratchb16/gavazzi/Images_HAGN/tmp"  # "/data61/laigle/Mock_lensing"
#    const DEFDIR="/scratchb16/gavazzi/Deflections_HAGN"
    const DEFDIR="/data62/gouin/DEFLECTION_PLANES"
    #const IMADIRS=[IMADIR]
    #const IMADIRS="/scratchb0".*string.(4:6).*"/laigle/Mock_lensing"
#    const IMADIRS=vcat( "/scratchb0".*string.(4:5).*"/laigle/Mock_lensing","/data78/laigle_scratchb06/Mock_lensing" )
    const IMADIRS=vcat( "/scratchb0".*string.([2,4,5,8]).*"/laigle/Mock_lensing",
                        "/scratchb21/laigle/Mock_lensing" )
    ind2slabs(id::Integer) = 20554 .- 40*[id-1,id]
    slab2ind(slid::Integer) = div(20554-slid,40)
    
    ########################################
    function ind2_CL_ima( id::Integer, band::AbstractString ; dust="dust", opa="" ) 
        for d in IMADIRS
            q=ind2slabs(id)
            t = joinpath( d, @sprintf( "Mock_sca_%s_%s%s_%05d_%05d_Init.dat", dust,
                                       opa, band, q[2], q[1] ) )
            isfile(t) && return t
            t=replace(t,r"_sca"=>"")
            isfile(t) && return t
        end
    end
    ########################################
    full_bin_name(id::Integer) = joinpath(DEFDIR,@sprintf( "test2_deflection_ipl_%04d_propage.bin", id ) )


    ########################################
"""
   read_bin( fn::AbstractString ; chunck::Int64=1_000_000, info::Bool=false )

Reads a deflection field map plus corresponding source redshift
# Arguments
- `fn`: input file  (Fortran binary)
- `chuncksize` : way file writing was done with Fortran records (shouldn't be changed) 
- `info` : if set to true, only reads header containing size and redshift

"""
    function read_bin( fn::AbstractString; chuncksize::Int64=1_000_000, info::Bool=false, reorder=true )
        ff = open( fn, "r" )
        f = FortranFile( ff )
        size = convert.( Int64, read( f, (Int32,2) ) )
        dum = read( f, (Float32,3) )
        zs = read( f, Float64 )
        if info
            close(ff)
            return (size, zs)
        end
        nt = prod(size)
        lck = min( nt, chuncksize )
        alpha = Array{Float32,3}(undef, size[1], size[2], 2 )
        i1=1
        dum = zeros( Float32, lck )
        while i1 <= nt
            i2 = min( i1+lck-1, nt )
            nv = i2-i1+1
            (nv<lck) && ( dum = zeros( Float32, nv ) )
            read( f, dum ) 
            @inbounds alpha[i1:i2] .= dum
            read( f, dum )
            @inbounds alpha[i1+nt:i2+nt] .= dum
            i1 += lck
        end
        close(ff)
        dum=nothing ; GC.gc()
        reorder || return (alpha, zs)
        alphat=permutedims(alpha,[3,1,2])
        alpha=nothing ; GC.gc()
        (alphat, zs)
    end
    ##############################################
    function read_bin_ima( fn::AbstractString )
        io=open(fn,"r")
        n=read(io,Int32)
        m=read(io,Int32)
        ima = Array{Float32}(undef, n, m)
        read!(io,ima)
        close(io)
        ima
    end
    #############################################
    wcs_header( z, dx, np ) = FITSHeader(  [ "REDSHIFT", "RADECSYS", "EQUINOX",
                                             "CTYPE1", "CUNIT1", "CRVAL1", "CRPIX1", "CDELT1", 
                                             "CTYPE2", "CUNIT2", "CRVAL2", "CRPIX2", "CDELT2" ] ,
                                           [ Float64(z), "ICRS", 2000.0,
                                             "RA---TAN", "deg", 0.0, (np[1]+1)/2, -dx/np[1],
                                             "DEC--TAN", "deg", 0.0, (np[2]+1)/2,  dx/np[2] ] ,
                                           [ "source plane redshift", "Astrometric system", "Mean equinox",
                                             "WCS projection type for this axis", "Axis unit", 
                                             "World coordinate on this axis", "Reference pixel on this axis",
                                             "Linear scaling",
                                             "WCS projection type for this axis", "Axis unit", 
                                             "World coordinate on this axis", "Reference pixel on this axis",
                                             "Linear scaling" ] 
                                           )

    ##################################################
    function ima_distort(id::Integer, band::AbstractString ; skip_deflection::Bool=false, dust="dust", opa="" )
        imaf=ind2_CL_ima(id,band,dust=dust,opa=opa)
        if !isfile(imaf)
            @warn "Cannot read in input image $imaf"
            return 0f0
        end
        @time ima0 = read_bin_ima(imaf)
        skip_deflection && return ima0
        fn=full_bin_name(id)
        if !isfile(fn) 
            @warn "Cannot read in deflection map $fn... returning unlensed image!"
            return ima0
        end
#        @time (defl_field,zz) = read_bin(fn,reorder=false)
        @time (defl_field,zz) = read_bin(fn,reorder=true)
        xmin, xmax, ymin, ymax = -dx/2.f0, dx/2.f0, -dx/2.f0, dx/2.f0
        xvals = range( xmin, xmax, length=size(ima0,1) )
        yvals = range( ymin, ymax, length=size(ima0,2) )
        @info "Worker $(myid()):  $(now()) For plane $id, read image in $imaf and deflection map in $fn"
        itp = LinearInterpolation( (collect(xvals),collect(yvals)), ima0, extrapolation_bc=0.f0 )
        ima = similar( ima0, axes(ima0) )
        Threads.@threads for j in axes(ima0,2)
            @inbounds for i in axes(ima0,2)
                ima[i,j] = itp( xvals[i] - defl_field[1,i,j],
                                yvals[j] - defl_field[2,i,j] )
#                ima[i,j] = itp( xvals[i] - defl_field[i,j,1],
#                                yvals[j] - defl_field[i,j,2] )
            end
        end
        ## never hurts to force release of memory
        ima0=itp=defl_field=nothing ; GC.gc()
        flush.((stdout,stderr))
        ima
    end
end #.... everywhere block

fils = ("u_LSST","VIS","Y_EUC","J_EUC","H_EUC")
#fils = split("u_LSST g_LSST r_LSST i_LSST z_LSST y_LSST VIS Y_EUC J_EUC H_EUC")

fil = length(ARGS)>0 ? ARGS[1] : fils[2]

dust="dust"  ## could be "nodust"
opa="opa_u_"  ## could be ""
#lensed=true  ## could be false to generate unlensed image
lensed=false  
#outdir=IMADIR
#outdir="/scratchb16/gavazzi/Images_HAGN"
outdir="/data41/gavazzi"

outfile = joinpath( outdir,
                    dust*"_"*opa*( lensed ? "lensed_" : "" )*"$fil.fits" )
#outfile=joinpath( outdir, "tmp.fits")


@info "Computations performed on $(gethostname()) with $(nworkers()) workers"*
       " and $(Threads.nthreads()) threads\nResult image for band $fil will be stored in $outfile"

# This is to get only existing CL images...
# Normally one would just write gids=2:499 or alike!
flist=String[]
for d in IMADIRS
    isdir(d) || continue
    for s in ("","sca_")
        pt=Regex("Mock_$(s)$(dust)_$(opa)$(fil)_\\d{5}_\\d{5}_Init.dat\$")
        append!(flist, filter( x->occursin(pt,x), readdir(d) ) )
    end
end
gids = slab2ind.( parse.(Int, replace.( flist, r".*_(\d{5})_(\d{5})_Init.dat"=>s"\1" ) ) )
sort!(gids) ; unique!(gids)
shuffle!(gids)
#gids=shuffle(gids)[1:6]
@info "Will combine $(length(gids)) planes: " *join(gids," ")
flush.((stdout,stderr))

if nprocs()>1
    gima = @distributed (+) for id in gids
        @info "$(now()) Reading slice $id"
        flush.((stdout,stderr))
        ima_distort( id, fil, skip_deflection=!lensed, dust=dust, opa=opa )
    end
else
    gima = ima_distort( first(gids), fil, skip_deflection=true, dust=dust, opa=opa ) .* 0.f0
    @showprogress 2.0 "Computing..." for id in gids
        gima .+= ima_distort( id, fil, skip_deflection=!lensed, dust=dust, opa=opa )
    end
end

@info "Done with computations... storing in $outfile"
io=FITS(outfile,"w")
write(io,gima,header=wcs_header(0.0f0,dx,size(gima)) )
close(io)
flush.((stdout,stderr))
