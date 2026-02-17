
using Distributed, Printf, SharedArrays,
    FortranFiles, FITSIO, Statistics,
    Interpolations, ProgressMeter, Plots, ForwardDiff
 
if occursin(r"amalgam", Base.Libc.gethostname())
    const dir0 = "/raid/gavazzi/ray-tracing-HAGN-bak"
    const tmpd = joinpath(dir0,"tmp_out")
    const dircat = dir0
else
    const dir0 = "/home/lmagri-stella/raytracing_test"
    const tmpd = "/net/GECO/nas12c/euclid/TestWorkFolder"
    const dircat = "/data75/laigle/CONE/Galaxies/catalogs"
end
const chunck_list = ["0-1","1-2","2-3","3-6"]
const oldgl = map( x->joinpath( dircat, "fits", "Galaxies_"*x*".fits"), chunck_list )
const gl = map( x->joinpath( dircat, "fits", "Galaxies_"*x*".d.fits"), chunck_list )
const hl = map( x->joinpath( dircat, "fits", "Haloes_"*x*".d.fits"), chunck_list )
const MatchMainF = map( x->joinpath( dircat, "MatchHaloGal", "Cat_"*x*"_Gal_MainHaloes.d.txt"), chunck_list )
const MatchSubF = map( x->joinpath( dircat, "MatchHaloGal", "Cat_"*x*"_Gal_SubHaloes.d.txt"), chunck_list )

zref = [ 0.2, 0.5, 0.8, 1.2, 2.0, 3.5, 5.0 ]

##################################################
"""
   average_law( x::Vector, y::Vector;  nbins=30, xrange=nothing, corrected=false )

Equivalent of average_loi old IDL funcion
# Arguments
- `x`: independent input array
- `y`: dependent input array
- `nbins` : number of bins
- `xrange` : min,max range over which `x` values are considered
- `corrected` : if set to true, unbiased variances are computed

"""
function average_law( x::Vector, y::Vector; nbins=30, xrange=nothing, corrected=false )

    ( xrange == nothing ) &&  ( xrange=extrema(x) )

    xbound = linspace(xrange[1],xrange[2],nbins+1)
    dx = xbound[2] - xbound[1]

    ym = zeros( eltype(x[1]), nbins)
    nc = zeros( Int64, nbins)
    ymed = zeros( eltype(x[1]), nbins)
    ydev = zeros( eltype(x[1]), nbins)
    for j=1:nbins
        ok = findall( xbound[j] .<= x .< xbound[j+1] )
        nc[j] =length(ok)
        if nc[j]>0
            ( ym[j], ydev[j] ) = StatsBase.mean_and_std( y[ok], corrected=corrected )
            ymed[j] = median( y[ok] )
        else
            ym[j]=NaN
            ymed[j]=NaN
            ydev[j]=NaN
        end
    end
    xbound = xbound[1:end-1]+dx/2

    return (xbound,ym,ymed,ydev,nc)
end


##############################################
"""
   read_bin( fn::AbstractString ; chunck::Int64=1_000_000, info::Bool=false )

Reads a binary deflection field map plus corresponding source redshift
# Arguments
- `fn`: File name where to write the model
- `chuncksize` : way file writing was done with Fortran records (shouldn't be changed) 
- `info` : if set to true, only reads header containing size and redshift

"""
function read_bin( fn::AbstractString; chuncksize::Int64=1_000_000, info::Bool=false )
    ff = try open( fn, "r" )
    catch
        @warn  "Cannot read in $fn"
        return ( nothing , nothing  )
    end
    f = FortranFile( ff )
    size = read( f, (Int32,2) ) ; size = convert.(Int64,size)
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
    return (alpha, zs)
end

##############################################
function read_bin_ima( fn::AbstractString )
    io=open(fn,"r")
    n=read(io,Int32)
    m=read(io,Int32)
    data = Array{Float32}(undef, n, m)
    read!(io,data)
    close(io)
    data
end


##############################################
"""
   deriv!(  x::Vector{Float32} , res::Vector{Float32} )

spatial derivatives (finite-differences) of a regularly-sampled array 
# Arguments
- `x`: input array
- `res` : resulting derivative (should be preallocated)

"""
##############################################
"""
function deriv!( x::Vector{Float32} , res::Vector{Float32} )
    n::Int64=length(x)
    res[1] = x[2]-x[1]
    @inbounds @simd for i=2:n-1
        res[i] = 0.5f0*( x[i+1] - x[i-1] )
    end
    res[n] = x[n]-x[n-1]
    nothing
end
"""
##############################################

function deriv!(x::Vector{Float32} , res::Vector{Float32})
	n::Int64=length(x)
	res[1]=-3*x[1]+4*x[2]-x[3]
	res[2]=-3*x[2]+4*x[3]-x[4]
	@inbounds @simd for i=3:n-2
		res[i]=-x[i-2]+8*x[i-1]-8*x[i+1]+x[i+2]
	end
	res[n-1]= 3*x[n-1]-4*x[n-2]+x[n-3]
	res[n]= 3*x[n]-4*x[n-1]+x[n-2]
	res = -res/12
	nothing
end
	

##############################################
"""

	Good one to use should be alpha2jac

"""
function alpha2jac_old( alpha::Array{Float32,3} ; scale::Float32=1.0f0 )
    nn = size(alpha)
    (nn1::Int64,nn2::Int64) = (nn[1],nn[2])
    (conv1,conv2) = (nn1/scale, nn2/scale)
    jac = Array{Float32,3}(undef,nn1,nn2,4)    
    d = Array{Float32,1}(undef, nn1)
    dd = Array{Float32,1}(undef, nn1)
    @inbounds for  j=1:nn2
        @inbounds @simd for  i=1:nn1
           d[i] = alpha[i,j,1] * conv1
        end
        deriv!( d, dd )
        @inbounds @simd for  i=1:nn1
            jac[i,j,1] = dd[i]
        end
    end
    @inbounds for  j=1:nn2
        @inbounds @simd for  i=1:nn1
           d[i] = alpha[i,j,2] * conv1
        end
        deriv!( d, dd )
        @inbounds @simd for  i=1:nn1
            jac[i,j,2] = dd[i]
        end
    end
    d = Array{Float32,1}(undef, nn2)
    dd = Array{Float32,1}(undef, nn2)
    @inbounds for i=1:nn1
        @inbounds @simd for  j=1:nn2
            d[j] = alpha[i,j,1] * conv2
        end
        deriv!( d, dd )
        @inbounds @simd for j=1:nn2
            jac[i,j,3] = dd[j]
        end
    end
    @inbounds for i=1:nn[1]
        @inbounds @simd for  j=1:nn2
            d[j] = alpha[i,j,2] * conv2
        end
        deriv!( d, dd )
        @inbounds @simd for j=1:nn2
            jac[i,j,4] = dd[j]
        end
    end 
    return jac
end

#############################################
# CRASH ZONE
#############################################

function alpha2jac( alpha::Array{Float32,3} ; scale::Float32=1.0f0 )
    nn = size(alpha)
    (nn1::Int64,nn2::Int64) = (nn[1],nn[2])
    (conv1,conv2) = (nn1/scale, nn2/scale)
    jac = Array{Float32,3}(undef,nn1,nn2,12)    
    d = Array{Float32,1}(undef, nn1)
    dd = Array{Float32,1}(undef, nn1)
    ddd = Array{Float32,1}(undef, nn1)

## Computing d1a1 and d1d1a1
    @inbounds for  j=1:nn2
        @inbounds @simd for  i=1:nn1
           d[i] = alpha[i,j,1] * conv1
        end
        deriv!( d, dd )
        deriv!( dd, ddd )
        @inbounds @simd for  i=1:nn1
            jac[i,j,1] = dd[i]
            jac[i,j,5] = ddd[i]
        end
    end

## Computing d1a2 and d1d1a2
    @inbounds for  j=1:nn2
        @inbounds @simd for  i=1:nn1
           d[i] = alpha[i,j,2] * conv1
        end
        deriv!( d, dd )
        deriv!( dd, ddd )
        @inbounds @simd for  i=1:nn1
            jac[i,j,2] = dd[i]
            jac[i,j,6] = ddd[i]
        end
    end

    d = Array{Float32,1}(undef, nn2)
    dd = Array{Float32,1}(undef, nn2)
    ddd = Array{Float32,1}(undef, nn2)

## Computing d2a1 and d2d2a1
    @inbounds for i=1:nn1
        @inbounds @simd for  j=1:nn2
            d[j] = alpha[i,j,1] * conv2
        end
        deriv!( d, dd )
        deriv!( dd, ddd )
        @inbounds @simd for j=1:nn2
            jac[i,j,3] = dd[j]
            jac[i,j,7] = ddd[j]
        end
    end

## Computing d2a2 and d2d2a2
    @inbounds for i=1:nn1
        @inbounds @simd for  j=1:nn2
            d[j] = alpha[i,j,2] * conv2
        end
        deriv!( d, dd )
        deriv!( dd, ddd )
        @inbounds @simd for j=1:nn2
            jac[i,j,4] = dd[j]
            jac[i,j,8] = ddd[j]
        end
    end

## Computing mixed derivatives

## d2d1a1
    @inbounds for i=1:nn1
        @inbounds @simd for  j=1:nn2
            d[j] = jac[i,j,1] /conv1
        end
        deriv!( d, dd )
        @inbounds @simd for j=1:nn2
            jac[i,j,9] = dd[j]
        end
    end
## d2d1a2
    @inbounds for i=1:nn1
        @inbounds @simd for  j=1:nn2
            d[j] = jac[i,j,2] /conv1
        end
        deriv!( d, dd )
        @inbounds @simd for j=1:nn2
            jac[i,j,10] = dd[j]
        end
    end
## d1d2a1
    @inbounds for j=1:nn2
        @inbounds @simd for i=1:nn1
            d[i] = jac[i,j,3] /conv2
        end
        deriv!( d, dd )
        @inbounds @simd for  i=1:nn1
            jac[i,j,11] = dd[j]
        end
    end
## d1d2a2
    @inbounds for j=1:nn2
        @inbounds @simd for i=1:nn1
            d[i] = jac[i,j,4] /conv2
        end
        deriv!( d, dd )
        @inbounds @simd for  i=1:nn1
            jac[i,j,12] = dd[j]
        end
    end
    return jac
end

##############################################
# alpha2jac test, considering commuting derivatives and no rotational terms (d1gamma2 = d2gamma1), validity to discuss
##############################################
function alpha2jac_norot_commder( alpha::Array{Float32,3} ; scale::Float32=1.0f0 )
    nn = size(alpha)
    (nn1::Int64,nn2::Int64) = (nn[1],nn[2])
    (conv1,conv2) = (nn1/scale, nn2/scale)
    jac = Array{Float32,3}(undef,nn1,nn2,8)    
    d = Array{Float32,1}(undef, nn1)
    dd = Array{Float32,1}(undef, nn1)
    ddd = Array{Float32,1}(undef, nn1)

    @inbounds for  j=1:nn2

        @inbounds @simd for  i=1:nn1
           d[i] = alpha[i,j,1] * conv1
        end

        deriv!( d, dd )
	deriv!( dd, ddd )

        @inbounds @simd for  i=1:nn1
            jac[i,j,1] = dd[i]
	    jac[i,j,5] = ddd[i]
        end

    end

    @inbounds for  j=1:nn2
        @inbounds @simd for  i=1:nn1
           d[i] = alpha[i,j,2] * conv1
        end
        deriv!( d, dd )
	deriv!( dd, ddd )
        @inbounds @simd for  i=1:nn1
            jac[i,j,2] = dd[i]
	    jac[i,j,6] = ddd[i]
        end
    end

    d = Array{Float32,1}(undef, nn2)
    dd = Array{Float32,1}(undef, nn2)
    ddd = Array{Float32,1}(undef, nn2)

    @inbounds for i=1:nn1
        @inbounds @simd for  j=1:nn2
            d[j] = alpha[i,j,1] * conv2
        end
        deriv!( d, dd )
	deriv!( dd, ddd )
        @inbounds @simd for j=1:nn2
            jac[i,j,3] = dd[j]
	    jac[i,j,7] = ddd[j]
        end
    end
    @inbounds for i=1:nn[1]
        @inbounds @simd for  j=1:nn2
            d[j] = alpha[i,j,2] * conv2
        end
        deriv!( d, dd )
	deriv!( dd, ddd )
        @inbounds @simd for j=1:nn2
            jac[i,j,4] = dd[j]
	    jac[i,j,8] = ddd[j]
        end
    end 
    return jac
end

##############################################
# Jacobian considering all terms, including rotationnal ones and the possibly non-commuting derivatives
##############################################

function alpha2jac_flexion_w_cross_terms( alpha::Array{Float32,3} ; scale::Float32=1.0f0 )
    nn = size(alpha)
    (nn1::Int64,nn2::Int64) = (nn[1],nn[2])
    (conv1,conv2) = (nn1/scale, nn2/scale)
    jac = Array{Float32,3}(undef,nn1,nn2,12)    
    d = Array{Float32,1}(undef, nn1)
    dd = Array{Float32,1}(undef, nn1)
    ddd = Array{Float32,1}(undef, nn1)

## This loop computes d1a1 and d1d1a1
    @inbounds for  j=1:nn2
        @inbounds @simd for  i=1:nn1
           d[i] = alpha[i,j,1] * conv1
        end
        deriv!( d, dd )
	deriv!( dd, ddd )
        @inbounds @simd for  i=1:nn1
            jac[i,j,1] = dd[i]
	    jac[i,j,5] = ddd[i]
        end
    end

## This loop computes d1a2 and d1d1a2
    @inbounds for  j=1:nn2
        @inbounds @simd for  i=1:nn1
           d[i] = alpha[i,j,2] * conv1
        end
        deriv!( d, dd )
	deriv!( dd, ddd )
        @inbounds @simd for  i=1:nn1
            jac[i,j,2] = dd[i]
	    jac[i,j,6] = ddd[i]
        end
    end

    d = Array{Float32,1}(undef, nn2)
    dd = Array{Float32,1}(undef, nn2)
    ddd = Array{Float32,1}(undef, nn2)

## This loop computes d2a1 and d2d2a1
    @inbounds for i=1:nn1
        @inbounds @simd for  j=1:nn2
            d[j] = alpha[i,j,1] * conv2
        end
        deriv!( d, dd )
	deriv!( dd, ddd )
        @inbounds @simd for j=1:nn2
            jac[i,j,3] = dd[j]
	    jac[i,j,7] = ddd[j]
        end
    end

## This loop computes d2a2 and d2d2a1
    @inbounds for i=1:nn[1]
        @inbounds @simd for  j=1:nn2
            d[j] = alpha[i,j,2] * conv2
        end
        deriv!( d, dd )
	deriv!( dd, ddd )
        @inbounds @simd for j=1:nn2
            jac[i,j,4] = dd[j]
	    jac[i,j,8] = ddd[j]
        end
    end

    d2a1 = Array{Float32,1}(undef,nn2)
    d2a2 = Array{Float32,1}(undef,nn2)
    d1d2a2 = Array{Float32,1}(undef,nn1)
    d1d2a1 = Array{Float32,1}(undef,nn1)
    
    d1a1 = Array{Float32,1}(undef,nn1)
    d1a2 = Array{Float32,1}(undef,nn1)
    d2d1a1 = Array{Float32,1}(undef,nn2)
    d2d1a2 = Array{Float32,1}(undef,nn2)
  
    d1a1grid = Array{Float32,2}(undef,nn1,nn2)
    d1a2grid = Array{Float32,2}(undef,nn1,nn2)
    d2a1grid = Array{Float32,2}(undef,nn1,nn2)
    d2a2grid = Array{Float32,2}(undef,nn1,nn2)

## These loops compute d1d2a1 and d1d2a2
    @inbounds for i=1:nn1
        @inbounds @simd for  j=1:nn2
            d2a1grid[i,j] = jac[i,j,3]
            d2a2grid[i,j] = jac[i,j,4]
        end
    end

    @inbounds for i=1:nn2
        @inbounds @simd for  j=1:nn1
            d2a1[j] = d2a1grid[j,i]
            d2a2[j] = d2a2grid[j,i]
        end
    

        deriv!( d2a1, d1d2a1 )
	deriv!( d2a2, d1d2a2 )
        @inbounds @simd for j=1:nn1
            jac[j,i,11] = d1d2a1[j]
	    jac[j,i,12] = d1d2a2[j]
        end
    end 
 
## These loops compute d2d1a1 and d2d1a2
    @inbounds for j=1:nn2
        @inbounds @simd for  i=1:nn1
            d1a1grid[i,j] = jac[i,j,1]
            d1a2grid[i,j] = jac[i,j,2]
            
        end
    end
      
    @inbounds for j=1:nn1
	@inbounds @simd for i=1:nn2
	    d1a1[i] = d1a1grid[j,i]
	    d1a2[i] = d1a1grid[j,i]
	end

        deriv!( d1a1, d2d1a1 )
	deriv!( d1a2, d2d1a2 )
        @inbounds @simd for i=1:nn2
            jac[j,i,9] = d2d1a1[i]
	    jac[j,i,10] = d2d1a2[i]
        end
    end    

    return jac
end
##############################################
# Another approach, taking first derivatives computed with alpha2jac and returning an array completed w/ next order derivatives
##############################################
function jac2der( jac::Array{Float32,3} )
	nn = size(jac)
    	(nn1::Int64,nn2::Int64) = (nn[1],nn[2])
	der = Array{Float32,3}(undef,nn1,nn2,8)
	d = Array{Float32,1}(undef, nn1)
	dd = Array{Float32,1}(undef, nn1)
	jacplusder = Array{Float32,3}(undef,nn1,nn2,12)
	jacplusder[:,:,1] = jac[:,:,1]
	jacplusder[:,:,2] = jac[:,:,2]
	jacplusder[:,:,3] = jac[:,:,3]
	jacplusder[:,:,4] = jac[:,:,4]

	@inbounds for j=1:nn2
		@inbounds @simd for i=1:nn1
			d[i] = jac[i,j,1] #selecting d1a1
		end

		deriv!(d,dd) #computing d1d1a1
		
		@inbounds @simd for i=1:nn1
			jacplusder[i,j,5] = d[i] #storing d1d1a1
		end
	end

	@inbounds for j=1:nn2
		@inbounds @simd for i=1:nn1
			d[i] = jac[i,j,2] #selecting d1a2
		end

		deriv!(d,dd) #computing d1d1a2
		
		@inbounds @simd for i=1:nn1
			jacplusder[i,j,6] = d[i] #storing d1d1a2
		end
	end

	@inbounds for j=1:nn2
		@inbounds @simd for i=1:nn1
			d[i] = jac[i,j,3] #selecting d2a1
		end

		deriv!(d,dd) #computing d1d2a1
		
		@inbounds @simd for i=1:nn1
			jacplusder[i,j,11] = d[i] #storing d1d2a1
		end
	end

	@inbounds for j=1:nn2
		@inbounds @simd for i=1:nn1
			d[i] = jac[i,j,4] #selecting d2a2
		end

		deriv!(d,dd) #computing d1d2a2
		
		@inbounds @simd for i=1:nn1
			jacplusder[i,j,12] = d[i] #storing d1d2a2
		end
	end

	d = Array{Float32,1}(undef, nn2)
	dd = Array{Float32,1}(undef, nn2)

	@inbounds for i=1:nn1
		@inbounds @simd for j=1:nn2
			d[j] = jac[i,j,3] #selecting d2a1
		end

		deriv!(d,dd) #computing d2d2a1
		
		@inbounds @simd for j=1:nn2
			jacplusder[i,j,7] = d[j] #storing d2d2a1
		end
	end

	@inbounds for i=1:nn1
		@inbounds @simd for j=1:nn2
			d[j] = jac[i,j,4] #selecting d2a2
		end

		deriv!(d,dd) #computing d2d2a2
		
		@inbounds @simd for j=1:nn2
			jacplusder[i,j,8] = d[j] #storing d2d2a2
		end
	end

	@inbounds for i=1:nn1
		@inbounds @simd for j=1:nn2
			d[j] = jac[i,j,1] #selecting d1a1
		end

		deriv!(d,dd) #computing d2d1a1
		
		@inbounds @simd for j=1:nn2
			jacplusder[i,j,9] = d[j] #storing d2d1a1
		end
	end


	@inbounds for i=1:nn1
		@inbounds @simd for j=1:nn2
			d[j] = jac[i,j,2] #selecting d1a2
		end

		deriv!(d,dd) #computing d2d1a2
		
		@inbounds @simd for j=1:nn2
			jacplusder[i,j,10] = d[j] #storing d2d1a2
		end
	end
	
	return jacplusder	
end

##############################################
function jac2muinv( jac::Array{Float32,3} ) 
    nn = size(jac)
    res = Array{Float32,2}( undef, nn[1], nn[2] )
    @inbounds for j=1:nn[2]
        @inbounds @simd for i=1:nn[1]
            res[i,j] = (1.0f0-jac[i,j,1])*(1.0f0-jac[i,j,4])-jac[i,j,2]*jac[i,j,3]
        end
    end
    return res
end
##############################################
function jac2kappa( jac::Array{Float32,3} ) 
    nn = size(jac)
    res = Array{Float32,2}( undef, nn[1], nn[2] )
    @inbounds for j=1:nn[2]
        @inbounds @simd for i=1:nn[1]
            res[i,j] = 0.5f0*( jac[i,j,1] + jac[i,j,4] )
        end
    end
    return res
end
##############################################
function jac2gamma1( jac::Array{Float32,3} ) 
    nn = size(jac)
    res = Array{Float32,2}( undef, nn[1], nn[2] )
    @inbounds for j=1:nn[2]
        @inbounds @simd for i=1:nn[1]
            res[i,j] = 0.5f0*( jac[i,j,1] - jac[i,j,4] )
        end
    end
    return res
end
##############################################
function jac2gamma2( jac::Array{Float32,3} ) 
    nn = size(jac)
    res = Array{Float32,2}( undef, nn[1], nn[2] )
    @inbounds for j=1:nn[2]
        @inbounds @simd for i=1:nn[1]
            res[i,j] = 0.5f0*( jac[i,j,2] + jac[i,j,3] )
        end
    end
    return res
end
##############################################
function jac2F1( jac::Array{Float32,3} ) 
    nn = size(jac)
    res = Array{Float32,2}( undef, nn[1], nn[2] )
    @inbounds for j=1:nn[2]
        @inbounds @simd for i=1:nn[1]
            res[i,j] = 0.5f0*(jac[i,j,5]-jac[i,j,12]) + 0.5f0*(jac[i,j,7]+jac[i,j,10])
        end
    end
    return res
end
##############################################
function jac2F2( jac::Array{Float32,3} ) 
    nn = size(jac)
    res = Array{Float32,2}( undef, nn[1], nn[2] )
    @inbounds for j=1:nn[2]
        @inbounds @simd for i=1:nn[1]
            res[i,j] = 0.5f0*(jac[i,j,6]+jac[i,j,11]) - 0.5f0*(jac[i,j,9]-jac[i,j,8])
        end
    end
    return res
end
##############################################
function jac2G1( jac::Array{Float32,3} ) 
    nn = size(jac)
    res = Array{Float32,2}( undef, nn[1], nn[2] )
    @inbounds for j=1:nn[2]
        @inbounds @simd for i=1:nn[1]
            res[i,j] = 0.5f0*(jac[i,j,5]-jac[i,j,12])-0.5f0*(jac[i,j,10]+jac[i,j,7])
        end
    end
    return res
end
##############################################
function jac2G2( jac::Array{Float32,3} )
    nn = size(jac)
    res = Array{Float32,2}( undef, nn[1], nn[2] )
    @inbounds for j=1:nn[2]
    	@inbounds @simd for i=1:nn[1]
	    res[i,j] = 0.5f0*(jac[i,j,6]+jac[i,j,11]) + 0.5f0*(jac[i,j,9]-jac[i,j,8])
	end
    end
    return res
end
##############################################
function jac2rot( jac::Array{Float32,3} ) 
    nn = size(jac)
    res = Array{Float32,2}( undef, nn[1], nn[2] )
    @inbounds for j=1:nn[2]
        @inbounds @simd for i=1:nn[1]
            res[i,j] = 0.5f0*( jac[i,j,2] - jac[i,j,3] )
        end
    end
    return res
end
##############################################
function jac2crossderivativea1( jac::Array{Float32,3} )
	nn = size(jac)
	res = Array{Float32,2}(undef, nn[1], nn[2] )
	@inbounds for j=1:nn[2]
		@inbounds @simd for i=1:nn[1]
			res[i,j] = 0.5f0*(jac[i,j,11]-jac[i,j,9])
		end
	end
	return res
end
##############################################
function jac2crossderivativea2( jac::Array{Float32,3} )
	nn = size(jac)
	res = Array{Float32,2}(undef, nn[1], nn[2] )
	@inbounds for j=1:nn[2]
		@inbounds @simd for i=1:nn[1]
			res[i,j] = 0.5f0*(jac[i,j,10]-jac[i,j,12])
		end
	end
	return res
end
##############################################
function alpha_mu_2_mus( alpha::Array{Float32,3}, sc )
    nn = size(alpha)
    p0 = ([nn[1], nn[2]] .+ 1.0f0)/2.0f0
    res = zeros( Float32, nn[1], nn[2] )
    rebin = 3
    hrb = div(rebin,2)
    w = (1.0f0/rebin)^2
    ssc=[sc[1],sc[2]]/rebin
    @inbounds for j=1:nn[2]
        y0=(j-p0[2]-0.5f0)*sc[2]
        #@show j
        @inbounds for i=1:nn[1]
            x0 = (i-p0[1]-0.5f0)*sc[1]
            #@show i,j
            if i>1 && i<nn[1] && j>1 && j<nn[2]
                @inbounds for jj=1:hrb
                    y = y0 + jj*ssc[2]
                    q = Float32(jj)/Float32(rebin)
                    @inbounds for ii=1:hrb
                        x = x0 + ii*ssc[1]
                        p = Float32(ii)/Float32(rebin)
                        a1 = alpha[i-1,j-1,1]*(0.5f0-q)*(0.5f0-p) + alpha[i,j-1,1]*(0.5f0-q)*(0.5f0+p) + 
                            alpha[i-1,j,1]*(0.5f0+q)*(0.5f0-p) + alpha[i,j,1]*(0.5f0+q)*(0.5f0+p)
                        a2 = alpha[i-1,j-1,2]*(0.5f0-q)*(0.5f0-p) + alpha[i,j-1,2]*(0.5f0-q)*(0.5f0+p) + 
                            alpha[i-1,j,2]*(0.5f0+q)*(0.5f0-p) + alpha[i,j,2]*(0.5f0+q)*(0.5f0+p)
                        u = x - a1
                        v = y - a2
                        ku = min( max( Int64( ceil( u/sc[1] + p0[1] ) ), 1 ), nn[1] )
                        kv = min( max( Int64( ceil( v/sc[2] + p0[2] ) ), 1 ), nn[2] )
                        res[ku,kv] += w
                    end
                    @inbounds for ii=hrb+1:rebin
                        x = x0 + ii*ssc[1]
                        p = Float32(ii)/Float32(rebin)
                        a1 = alpha[i,j-1,1]*(0.5f0-q)*(1.5f0-p) + alpha[i+1,j-1,1]*(0.5f0-q)*(p-0.5f0) + 
                            alpha[i,j,1]*(0.5f0+q)*(1.5f0-p) + alpha[i+1,j,1]*(0.5f0+q)*(p-0.5f0)
                        a2 = alpha[i,j-1,2]*(0.5f0-q)*(1.5f0-p) + alpha[i+1,j-1,2]*(0.5f0-q)*(p-0.5f0) + 
                            alpha[i,j,2]*(0.5f0+q)*(1.5f0-p) + alpha[i+1,j,2]*(0.5f0+q)*(p-0.5f0)
                        u = x - a1
                        v = y - a2
                        ku = min( max( Int64( ceil( u/sc[1] + p0[1] ) ), 1 ), nn[1] )
                        kv = min( max( Int64( ceil( v/sc[2] + p0[2] ) ), 1 ), nn[2] )
                        #@show "in ",ii,jj,ku,kv
                        res[ku,kv] += w
                    end
                end
                @inbounds for jj=hrb+1:rebin
                    y = y0 + jj*ssc[2]
                    q = Float32(jj)/Float32(rebin)
                    @inbounds for ii=1:hrb
                        x = x0 + ii*ssc[1]
                        p = Float32(ii)/Float32(rebin)
                        a1 = alpha[i-1,j,1]*(1.5f0-q)*(0.5f0-p) + alpha[i,j,1]*(1.5f0-q)*(0.5f0+p) + 
                            alpha[i-1,j+1,1]*(q-0.5f0)*(0.5f0-p) + alpha[i,j+1,1]*(q-0.5f0)*(0.5f0+p)
                        a2 = alpha[i-1,j,2]*(1.5f0-q)*(0.5f0-p) + alpha[i,j,2]*(1.5f0-q)*(0.5f0+p) + 
                            alpha[i-1,j+1,2]*(q-0.5f0)*(0.5f0-p) + alpha[i,j+1,2]*(q-0.5f0)*(0.5f0+p)
                        u = x - a1
                        v = y - a2
                        ku = min( max( Int64( ceil( u/sc[1] + p0[1] ) ), 1 ), nn[1] )
                        kv = min( max( Int64( ceil( v/sc[2] + p0[2] ) ), 1 ), nn[2] )
                        #@show "in ",ii,jj,ku,kv
                        res[ku,kv] += w
                    end
                    @inbounds for ii=hrb+1:rebin
                        x = x0 + ii*ssc[1]
                        p = Float32(ii)/Float32(rebin)
                        a1 = alpha[i,j,1]*(1.5f0-q)*(1.5f0-p) + alpha[i+1,j,1]*(1.5f0-q)*(p-0.5f0) + 
                            alpha[i,j+1,1]*(q-0.5f0)*(1.5f0-p) + alpha[i+1,j+1,1]*(q-0.5f0)*(p-0.5f0)
                        a2 = alpha[i,j,2]*(1.5f0-q)*(1.5f0-p) + alpha[i+1,j,2]*(1.5f0-q)*(p-0.5f0) + 
                            alpha[i,j+1,2]*(q-0.5f0)*(1.5f0-p) + alpha[i+1,j+1,2]*(q-0.5f0)*(p-0.5f0)
                        u = x - a1
                        v = y - a2
                        ku = min( max( Int64( ceil( u/sc[1] + p0[1] ) ), 1 ), nn[1] )
                        kv = min( max( Int64( ceil( v/sc[2] + p0[2] ) ), 1 ), nn[2] )
                        #@show "in ",ii,jj,ku,kv
                        res[ku,kv] += w
                    end
                end
            else
                x=x0+0.5*sc[1]
                y=y0+0.5*sc[2]
                u = x - alpha[i,j,1]
                v = y - alpha[i,j,2]
                ku = Int64( ceil( u/sc[1] + p0[1] ) )
                kv = Int64( ceil( v/sc[2] + p0[2] ) )
                0<ku<=nn[1] && 0<kv<=nn[2]  && ( res[ku,kv] += 1 )
            end


        end
    end
    return res
end

#############################################
function get_id_z_mapping( ; show=false, 
                           sd=occursin(r"amalgam", Base.Libc.gethostname()) ? "." : dir0*"/small_angle_OBB" )
    res = Dict{Int64, Float64}()
    if occursin( r"small_angle_SPL", sd )
        re = r".*_(\d+)_propage\.bin"
    elseif occursin( r"large_angle_SPL", sd )
        re = r".*_(\d+)_propage\.bin"
    elseif occursin( r"small_angle_OBB", sd )
        re = r".*_(\d+)_propage\.bin"
    else
        re = r".*_(\d+)\.bin"
    end
    for (root,dirs,files) in walkdir( sd, follow_symlinks=true )
        for f in files[1:end]
            if occursin(r".*\.bin",f)
                (s,z)=read_bin( joinpath(root,f), info=true )
                ff = replace( f, re => s"\1" )
                show && ( @show f,ff,z )
                id = parse(Int64, ff)
                res[id] = z
            end
        end
    end
    kmp = sort(collect( keys(res)))
    return ( kmp, map( x-> res[x], kmp) )
end

#############################################
function z2id( zin )
    (id,z) = get_id_z_mapping( )
    itp = interpolate( (z,), id, Gridded(Constant()))
    return Int64.( itp(zin) )
end
#############################################
function id2z( id0 )
    (id,z) = get_id_z_mapping( )
    itp = interpolate( (id,), z, Gridded(Constant()))
    return ( itp(id0) )
end

#############################################
function wcs_header( z, dx, np )
    FITSHeader( 
        [ "REDSHIFT", "RADECSYS", "EQUINOX",
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
end


#############################################
function full_bin_name( id::Int64 ; mode::AbstractString="OBB", size::Number=1.0 )
    occursin( r"amalgam", Base.Libc.gethostname())  && 
        return joinpath(dir0,@sprintf("test2_deflection_ipl_%04d.bin", id))
    if mode == "OBB" 
        if size<1.1 
            return @sprintf( "/net/GECO/nas12c/euclid/HAGN_lightcone/test2_deflection_ipl_%04d_propage.bin", id)
            #return @sprintf( "/data62/gouin/DEFLECTION_PLANES/test2_deflection_ipl_%04d_propage.bin", id)
        else
            return @sprintf( "/net/GECO/nas12c/euclid/HAGN_lightcone/test2_deflection_ipl_%04d_propage.bin", id)
        end
    elseif mode=="SPL"
        if size<1.1 
            return @sprintf( "/scratchb03/gouin/DEFLECTION_PLANES_SPL2/propagation_alpha_05_pad1_%04d_propage.bin", id )
        else
            return @sprintf( "/data62/gouin/DEFLECTION_PLANES_SPL2/propagation_alpha_1125_pad1_%04d_propage.bin", id)
        end
    else
        error("wrong mode  input. Must be either OBB or SPL")
    end
end


##############################################
function find_img( aa::Array{Float32,3}, u::Vector{Float32}, lev::Int64, res::Vector{Tuple{Int64,Int64}}, 
                   ilow::Int64, ihigh::Int64, jlow::Int64, jhigh::Int64  ; tol=100 )

    cc = inbox( aa[ilow,jlow,:], aa[ihigh,jlow,:], aa[ihigh,jhigh,:] ,aa[ilow,jhigh,:], u )
    #    ( !cc && lev==1 ) &&  return  ### outside at first level... we give up
    cc || return  ### outside at first level... we give up
    
    dist = max( abs(ihigh-ilow), abs(jhigh-jlow) )
    #    @show  "passed", lev, ilow, ihigh, jlow, jhigh , u, dist, tol, dist < tol

    if (dist < tol) ## || ( !cc && lev>1 ) ## sufficiently close or apparently outside
        @inbounds for j=jlow:jhigh-1
            @inbounds for i=ilow:ihigh-1
                ( intriangle(  aa[i,j,:], aa[i+1,j,:] ,  aa[i,j+1,:] , u ) |
                  intriangle(  aa[i+1,j,:], aa[i+1,j+1,:] ,  aa[i,j+1,:] , u )  ) &&
                  push!( res, (i,j) ) 
            end
        end
    else
        imp=div(ihigh+ilow,2)
        jmp=div(jhigh+jlow,2)
        find_img( aa, u, lev+1, res, ilow, imp, jlow, jmp, tol=tol )
        find_img( aa, u, lev+1, res, ilow, imp, jmp+1, jhigh , tol=tol )
        find_img( aa, u, lev+1, res, imp+1, ihigh, jlow, jmp , tol=tol)
        find_img( aa, u, lev+1, res, imp+1, ihigh, jmp+1, jhigh, tol=tol )
    end
end

#####################################################
function intriangle( pa::Vector{Float32}, pb::Vector{Float32},
                     pc::Vector{Float32}, px::Vector{Float32} )
#    const tol=-3f-6
#    const tol=-3f-4
    p=px[1]-pa[1]
    q=px[2]-pa[2]
    r=pb[1]-pa[1]
    s=pb[2]-pa[2]
    u=pc[1]-pa[1]
    v=pc[2]-pa[2]
    a=r*v-u*s  ; a2=a*a
    b=p*v-q*u
    d=q*r-p*s
    e=a-b-d
#    return all( [ a2-b*b, a2-d*d, a2-e*e, a*b, a*d, a*e ] .>= tol*a2 )
    return all( [ a2-b*b, a2-d*d, a2-e*e, a*b, a*d, a*e ] .>= -3f-6 *a2 )
end 
#####################################################
function inbox( pa::Vector{Float32}, pb::Vector{Float32}, 
                pc::Vector{Float32}, pd::Vector{Float32}, 
                px::Vector{Float32} )
    xr = extrema( [pa[1] , pb[1], pc[1], pd[1] ] )
    yr = extrema( [pa[2] , pb[2], pc[2], pd[2] ] )
    return  (xr[1] <= px[1] <= xr[2]) & (yr[1] <= px[2] <= yr[2])
end 


#####################################################
function read_all_cat_old( ; rmax=0.51f0 )
    id=nothing 
    ra=nothing
    dec=nothing
    ms=nothing
    z=nothing
    offset=0 
    for (il,gf) in enumerate(oldgl)        
        t=FITS(gf)
        id0=Int64.(read(t[2],"ID"))
        mm=maximum(id0)
        id0 .+=  offset
        ra0=-read(t[2],"Ra")
        dec0=read(t[2],"Dec")
        ms0=read(t[2],"Mass")
        z0=read(t[2],"z")
        ### sort with increasing redshift!!
        sp=sortperm(z0)
        id0=id0[sp]
        ra0=ra0[sp]
        dec0=dec0[sp]
        ms0=ms0[sp]
        z0=z0[sp]
        ### chop off wide fov stuff unless rmax is increased
        sp1 = findall( (-rmax .< ra0 .<= rmax) .&  (-rmax .< dec0 .<= rmax)  )
        id0=id0[sp1]
        ra0=ra0[sp1]
        dec0=dec0[sp1]
        ms0=ms0[sp1]
        z0=z0[sp1]
        if il==1 
            id=id0
            ra=ra0
            dec=dec0
            ms=ms0
            z = z0
        else
            append!(id,id0)
            append!(ra,ra0)
            append!(dec,dec0)
            append!(ms,ms0)
            append!(z,z0)
        end
        close(t)
        offset += mm
    end
    return (id,ra,dec,z,ms)
end
#####################################################
### Gal_MS ID=9, Halo_Mh ID=10, Halo_Mvir ID=22 
### Gal_a ID=56 Gal_b ID=57 Gal_et ID=58
function read_new_cat( gf ; sort=:z, rmax=0.501, massid=9 ,shapeid=-1 )  
    t = FITS(gf)
    gal = read( t[1] )
    close(t)
    id0 = Int64.( gal[1,:] )
    ra0 = Float64.( -gal[2,:] )
    dec0 = Float64.( gal[3,:] )
    z0 = Float64.( gal[4,:] )
    ms0 = Float64.( gal[massid,:] )
    if shapeid>0 
        a_ell = Float64.( gal[shapeid,:] )
        b_ell = Float64.( gal[shapeid+1,:] )
        th_ell = Float64.( gal[shapeid+2,:] )
        e_ell = ( a_ell .- b_ell ) ./ ( a_ell .+ b_ell )
        e1 = e_ell .* cos.( 2.0*th_ell )
        e2 = e_ell .* sin.( 2.0*th_ell )
    end

    if sort == :z    ### sort with increasing redshift!!
        sp = sortperm(z0)
    elseif sort == :id
        sp = sortperm(id0)
    else
        sp = 1:length(id0)
    end
    id0=id0[sp]
    ra0=ra0[sp]
    dec0=dec0[sp]
    ms0=ms0[sp]
    z0=z0[sp]
    if shapeid>0 
        e1= e1[sp]
        e2= e2[sp]
    end
    ### chop off wide fov stuff unless rmax is increased
    sp1 = findall( (-rmax .< ra0 .<= rmax) .&  (-rmax .< dec0 .<= rmax)  )
    id0=id0[sp1]
    ra0=ra0[sp1]
    dec0=dec0[sp1]
    ms0=ms0[sp1]
    z0=z0[sp1]
    if shapeid>0 
        e1= e1[sp1]
        e2= e2[sp1]
        return (id0,ra0,dec0,z0,ms0,e1,e2)
    else
        return (id0,ra0,dec0,z0,ms0)
    end
end

##############################################
##############################################
##############################################


#####################################################
function dump_imu_map( z ; dx=1.0 , mode="OBB" )
    ii = z2id( z ) 
    sii = @sprintf("%04d",ii) 
    fn = full_bin_name( ii, mode=mode )
    @show sii,z,fn
    (a,zz) = read_bin( fn )
    np = size(a)
    @show zz, z
    
    jac = alpha2jac( a, scale=Float32(dx) )
    @time imu = jac2muinv( jac )
    f = FITS( joinpath(tmpd,"imu-$mode-$sii.fits"), "w");  
    write( f, imu, header=wcs_header( zz, dx, np ) )
    close(f)
    gc(true)
end
#####################################################
function dump_kappa_map( z ; dx=1.0 , mode="OBB" , margin=0 )
    ii = z2id( z ) 
    sii = @sprintf("%04d",ii)
    fn = full_bin_name( ii, mode=mode )
    @show sii,z,fn
    (a,zz) = read_bin( fn )
    np = size(a)
    jac = alpha2jac( a, scale=Float32(dx) )
    @time kappa = jac2kappa( jac )
    f = FITS( joinpath(tmpd,"kappa-$mode-$sii.fits"), "w");  
    write( f, kappa, header=wcs_header( zz, dx, np ) )
    close(f)
    cs=(np[1]-margin*2)*(np[2]-margin*2)

    tmpk = reshape(kappa[margin+1:np[1]-margin, margin+1:np[2]-margin], cs)
    
    mu = mean( tmpk )
    med = median( tmpk )
    dev = std( tmpk )
    @printf("CVG   mean=%.4f median=%.4f  rms=%.4f   MODE=%s\n",mu,med,dev,mode)
    tit = @sprintf(" source redshift=%.3f",zz)
    @show tit, size(tmpk), cs

    if margin>0 
        @time histogram( tmpk,
                   title=tit, label="PDF", xlabel="convergence",  yscale=:log10, 
                   ylabel="Probability density", normalize=true , 
                   xlims=(-0.3,1.2), ylims=(1e-3,30) )
        vline!( [mu], label=@sprintf("mean = %.3f",mu) )
        vline!( [med], label=@sprintf("median = %.3f",med)  )
        vline!( mu+dev*[-1,1] , label=@sprintf("rms = %.3f",dev) )
        @time savefig( joinpath(tmpd,"kappa_PDF_$mode-$sii.pdf") )
    end
    gc(true)
end
##################################################
function dump_caustic_map( ii, z ; do_img_tau=false, x_section_pdf=false, 
                                       dx=1.0, mode="OBB", margin=0, mu_min=10.0 )

#    ii = z2id( z ) 
    sii = @sprintf("%04d",ii)
    fn = full_bin_name( ii, mode=mode, size=dx )

    @show ii, z, mode, fn

    (a,zz) = read_bin( fn )
    np = size(a)
    ps = dx ./ collect(np)

    nnp = collect(np) .- 2.0*margin
    cs = nnp[1]*nnp[2]

    @time cum = alpha_mu_2_mus( a, ps )
#    cum = zeros( Float32, np[1], np[2] )
    scum = view( cum, margin+1:np[1]-margin, margin+1:np[2]-margin )
#    scum[:,:] = randn(nnp[1],nnp[2]) 

    cfits = joinpath( tmpd, "caustics-$mode-$sii.fits" )
    f = FITS( cfits, "w" );  
    write( f, cum[ margin+1:np[1]-margin, margin+1:np[2]-margin], 
           header=wcs_header( zz, dx*(1-2.0*margin/np[1]), nnp ) )
    close( f )

    if do_img_tau
        jac = alpha2jac( a, scale=Float32(dx) )
        @time imu = jac2muinv( jac )
        jac = nothing     ## free memory
        imu = abs.(imu)
        simu = view( imu, margin+1:np[1]-margin, margin+1:np[2]-margin )
        sumsimu = sum(simu)
    end

    a = nothing    ## free memory

    nmu = 15
    imuthresh = 10 .^range( -log10(30.0) , stop=-log10(2.0), length=nmu )
    tau_s_true = zeros(nmu)

    ( tau_i, tau_s ) = ( nothing, nothing )
    do_img_tau && ( ( tau_i, tau_s ) = ( zeros(nmu), zeros(nmu) ) )
        tau_s = zeros(nmu)

    #### Dumps tau results
    f = open( joinpath( tmpd, "depth-$mode-$sii.txt" ), "w"  )
    print(f,"#!       mu      tau_s_true")
    do_img_tau &&  print(f,"    tau_i   tau_s")
    print(f,"\n")
    for i=1:nmu
        if do_img_tau
            ok = simu .<=  imuthresh[i]
            tau_i[i] = count( ok ) / cs
            tau_s[i] = sum( simu[ok] ) / sumsimu
        end
        ok = scum .>  1.0/imuthresh[i]
        tau_s_true[i] = count( ok ) / cs
        @printf(f,"   %.5f    %.5f", 1.0/imuthresh[i], tau_s_true[i] )
        do_img_tau && ( @printf(f,"   %.5f    %.5f", tau_i[i], tau_s[i] ) )
        print(f,"\n")
    end
    close(f)
    simu=nothing
    do_img_tau && (imu=nothing)
    scum = nothing
    cum = nothing
    gc(true)
    
    if x_section_pdf
        ccat = joinpath( tmpd, "caustics-$mode-$sii.txt" )
        scom = occursin( r"amalgam", Base.Libc.gethostname()) ? "sex"  :
            "/softs/sextractor/2.19.5-CentOS6/bin/sex"
        run(`$scom -c $dir0/sex.cfg  $cfits  -CATALOG_NAME  $ccat -DETECT_THRESH  $mu_min`)
    end
    
    return ( 1.0/imuthresh, tau_s_true, nothing, nothing )

end
#####################################################
function slow_src_img_catalog_mapping()
    t=FITS(gl[1])
    id=read(t[2],"ID")
    ra=-read(t[2],"Ra")
    dec=read(t[2],"Dec")
    ms=read(t[2],"Mass")
    z=read(t[2],"z")
    close(t)

    mrad=0.1/3600.
    dx=1.0
    fi = open( joinpath( tmpd, "cati.csv" ), "w"  )
    fs = open( joinpath( tmpd, "cats.csv" ), "w"  )
    for iz=1:length(zall)
        zu=zall[iz]
        zl= iz==1 ? 0.0 : zall[iz-1]
        zu > 1.0 && continue
        zu < 0.8 && continue
        @show zall[iz],idall[iz]
        ok = findall( zl .< z .<= zu )
        @show length(ok)
        println( fi, "######## $zl  $zu ", length(ok) )
        println( fs, "######## $zl  $zu ", length(ok) )
        for mid in ok
            @printf( fs, "%d  %f %f %f %f %f %f %f %f\n",
                     id[mid],ra[mid],dec[mid],ra[mid],dec[mid],
                     ra[mid],dec[mid],z[mid],ms[mid] )
        end
        
        @time tr = NearestNeighbors.KDTree( [ra[ok] dec[ok]]', Euclidean() ; )
        @show tr
        @time  (aa,zz) = read_bin( full_bin_name( idall[iz], mode="OBB" ) )
        nn=size(aa)
        cd1 = -dx/nn[1]
        cr1 = (nn[1]+1)/2
        cd2 = dx/nn[2]
        cr2 = (nn[2]+1)/2
        d=1.0/nn[1]
        @show nn,d,zz

        @inbounds for j=1:nn[2]
#        @showprogress 2 "Loop over grid rows..." for j=1:nn[2]
            y = cd2 *(j-cr2)
            #        @printf(STDERR,"\r%d/%d",j,nn[2])
            @inbounds for i=1:nn[1]
                x = cd1 *(i-cr1)
                u =  x + aa[i,j,1]
                v =  y - aa[i,j,2]
                res = inrange( tr, [u,v], mrad, true )
                if length(res)>0
                    mid=ok[res[1]]
                    @printf( fi, "%d  %f %f %f %f %f %f %f %f\n",
                             id[mid],ra[mid],dec[mid],x,y,
                             u,v,z[mid],ms[mid] )
                end
            end
        end
        print(STDERR,"\n")
        a=nothing
        tr=nothing
    end
    close(fi)
    close(fs)
end
#####################################################
function get_bbox( x, y, z, cr1, cr2, cd1, cd2, nn )
    res = Vector{Int64}(undef,4)
    W = Int64( floor( 400 + (z-0.05)*(1200-400)/(5-0.05) ) )
    i0 = Int64( floor( x/cd1 + cr1 ) )
    j0 = Int64( floor( y/cd2 + cr2 ) )
    res[1] = max(1,i0-W)
    res[2] = min(nn[1],i0+W)
    res[3] = max(1,j0-W)
    res[4] = min(nn[2],j0+W)
    return res
end
#####################################################
in_margin(x,y,tm) = ( abs(x)<tm ) & ( abs(y)<tm )
#####################################################
function fast_src_img_catalog_mapping( iz::Int64, croot::AbstractString , id, ra, dec, zcat, ms ;
                                       mode="OBB", dx=1.0f0, shuffle=:none, born=false )
    zlow = iz==1 ? 0.f0 : zall[iz-1]  ;   zhigh = zall[iz]
#    zlow = zall[iz] ;  zhigh= iz==endof(zall) ? 10. : zall[iz+1]    
    ok = findall( zlow .<= zcat .< zhigh )
    two_margin = dx - 2.f0 * 0.02f0  ## 1.2' thick
    fn = full_bin_name( idall[iz], mode=mode, size=dx )
    cn = @sprintf("%s_slice_%04d.cat",croot,idall[iz] ) 
    ( Float32(dx) > 1.1f0 ) && (cn = "wide_$cn")  ;  (mode=="OBB") || (cn = "SPL_$cn")
    @show length(ok), zlow, zhigh, fn, cn
    ex = false
    f = nothing
    try 
        f = open( joinpath( tmpd, cn ), "r"  )
        ex = true
        info("file already exists... skipping...")
        close(f)
    catch
    end
    ex && return 
    if length(ok)>0 
        @time (aa,zz) = read_bin( fn )
        @info "done reading deflection map!" 
        nn=size(aa)
        @time jac=alpha2jac( aa, scale=Float32(dx) )
        @info "done computing jacobian matrix map" 
        cd1 = -dx/nn[1]
        cr1 = (nn[1]+1)/2f0
        cd2 = dx/nn[2]
        cr2 = (nn[2]+1)/2f0
        @time begin
            x=0.f0 ; y=0.f0
            @inbounds @simd for i=1:nn[1]
                x = cd1 *(i-cr1)
                aa[i,:,1] .+= x    ### this is now u (negated)
            end
            @inbounds @simd for j=1:nn[2]
                y = cd2 *(j-cr2)
                aa[:,j,2] .= y .- aa[:,j,2]   ### this is now v
            end
        end
        @info "done turning alpha map into source plane positions"
    end
    N_tol_Max=4 ; tol0=40.0f0
    f = open( joinpath( tmpd, cn ), "w"  )
    print(f, "#!     id     ra_img    dec_img  multi  ra_src    dec_src      z       ms      a11      a12      a21     a22\n")
    for i in ok
        isolve=1
        du = [ 0.0f0, 0.0f0 ]
        u = (shuffle==:none) ? 
            convert.( Float32, [ra[i], dec[i] ] ) : 
            dx*(rand(Float32,2)-0.5f0)
        j = [ floor(Int64,u[1]/cd1+cr1), floor(Int64,u[2]/cd2+cr2) ]
        j[1] = max( 1, min(j[1],nn[1]) )  ;  j[2] = max( 1, min(j[2],nn[2]) )

        res = born ? [ (j[1],j[2]) ] : Vector{Tuple{Int64,Int64}}()
        while length(res)<1  && in_margin( u[1], u[2], two_margin ) && isolve<N_tol_Max 
            tol = tol0 * (1.0f0 + 4.f0 *(isolve-1))
            bbox = get_bbox( u[1]+du[1], u[2]+du[2], zall[iz], cr1, cr2, cd1, cd2, nn )
            find_img( aa, u, 1, res, bbox[1], bbox[2], bbox[3], bbox[4], tol=tol )
            du = (rand(Float32,2) .- 0.5f0) .* tol0 .* [cd1,cd2]
            isolve += 1
        end
        for (im,j) in enumerate(res)
            @printf(f, "   %d  %f  %f  %d  %f  %f  %f  %f   %f  %f  %f  %f\n",
                    id[i], cd1*(j[1]-cr1), cd2*(j[2]-cr2), im, u[1], u[2],
                    zcat[i], ms[i] ,
                    jac[j[1],j[2],1], jac[j[1],j[2],2],
                    jac[j[1],j[2],3], jac[j[1],j[2],4] )
        end
        (length(res)<1 ) && 
            @printf(f, "    %d   -----   -----   %d  %f  %f  %f  %f   -----   -----   -----   -----\n",
                    id[i], 0, u[1], u[2], zcat[i], ms[i] )
    end
    close(f)  ; @info "done storing new catalog"
    aa=nothing
    jac=nothing
    ok=nothing
    GC.gc(true)
end



##############################################
##############################################
##############################################
#############################################
function read_one_cat( fn::AbstractString ; born=false )
    a = try readdlm(fn, skipstart=1)
    catch
        return Tuple( fill(nothing,9) ) 
    end
    ID = convert.( Int64, a[:,1] )
    if count( ID .< 0 )>0
        @show fn, count( ID .< 0 ), size(a), typeof(a)
    end
    if born
        x = map( x -> isa(x,Float64) ? x : NaN , a[:,5] )
        y = map( x -> isa(x,Float64) ? x : NaN , a[:,6] )
    else
        x = map( x -> isa(x,Float64) ? x : NaN , a[:,2] )
        y = map( x -> isa(x,Float64) ? x : NaN , a[:,3] )
    end
    z = map( x -> isa(x,Float64) ? x : NaN , a[:,7] )
    ms = map( x -> isa(x,Float64) ? x : NaN , a[:,8] )
    a11 = map( x -> isa(x,Float64) ? x : NaN , a[:,9] )
    a12 = map( x -> isa(x,Float64) ? x : NaN , a[:,10] )
    a21 = map( x -> isa(x,Float64) ? x : NaN , a[:,11] )
    a22 = map( x -> isa(x,Float64) ? x : NaN , a[:,12] )
    imu = (1.0 .- a11) .* (1.0 .- a22) .-  a12 .* a21
    kappa  = 0.5 * ( a11 .+ a22 )
    gamma1 = 0.5 * ( a11 .- a22 )
    gamma2 = 0.5 * ( a21 .+ a12 )
    ok=findall( .!(isnan.(x)) )
    return ( ID[ok], x[ok], y[ok], z[ok], ms[ok], 1 ./ Array(imu')[ok], kappa[ok], gamma1[ok], gamma2[ok] )
end
#############################################
function read_all_cat( ; re::Regex=r"Galaxies_.*_slice_.*\\.cat" , born=false )
    (ID,x,y,z,ms,mu,k,g1,g2) = Tuple( fill(nothing,9) ) 
    for (root,dirs,files) in walkdir( tmpd ) 
        nf=0
        for f in files[1:end]
            if occursin(re,f)
@show "ok", f, joinpath(root,f)
                (IDt,xt,yt,zt,mst,mut,kt,g1t,g2t)= read_one_cat( joinpath(root,f) , born=born )
                @show length(xt)
                xt == nothing && continue
                nf += 1
                if nf == 1 
                    ID=IDt ; x=xt ; y=yt ; z=zt ; ms=mst ; mu=mut ; k=kt ; g1=g1t ; g2=g2t
                else
                    append!(ID,IDt)
                    append!(x,xt)
                    append!(y,yt)
                    append!(z,zt)
                    append!(ms,mst)
                    append!(mu,mut)
                    append!(k,kt)
                    append!(g1,g1t)
                    append!(g2,g2t)
                end
              #  @show f,nf,size(x)
            end
        end
    end
    return (ID,x,y,z,ms,mu,k,g1,g2)
end


#############################################
function build_cat( chk::Integer ; root="", rmax=0.501, born=false )

    ( (chk<1) || (chk>length(chunck_list)) ) &&  (error("wrong chunck id in build_cat"))

    ch=chunck_list[chk]

    (IDcat,xcat,ycat,zcat,mscat,mucat,kappacat,g1cat,g2cat) = 
        read_all_cat( re=Regex("^Galaxies_"*ch*"_slice_.*cat"), born=born )
#        read_all_cat( re=Regex(root*"Galaxies_"*ch*"_slice_.*\\.cat"), born=born )
    ng=length(IDcat)
    @show "lensed galaxies",ng, IDcat[ [1,end] ] , IDcat[end]-IDcat[1]+1

    (id0t,ra0t,dec0t,z0t,ms0t,e1t,e2t) = read_new_cat( gl[chk], sort=:none, rmax=100.0f0, shapeid=56 ) 
    nht=length(id0t) 
    @show "ref galaxies", nht, id0t[ [1,end] ] , extrema( id0t ) 

    (id0,ra0,dec0,z0,ms0) = read_new_cat( hl[chk], sort=:none, massid=10, rmax=100.0f0 ) 
    nh=length(id0) 
    @show "ref halos", nh, id0[ [1,end] ] ,  extrema( id0 ) 

    Mh_main = fill( NaN, ng )
    Mh_sub = fill( NaN, ng )
    es1_cat = fill( NaN, ng )
    es2_cat = fill( NaN, ng )
    e1_cat = fill( NaN, ng )
    e2_cat = fill( NaN, ng )

    sp1 = findall( (-rmax .< ra0t .<= rmax) .&  (-rmax .< dec0t .<= rmax)  )
    id0t = id0t[ sp1 ]
    if chk<4
        as=convert.(Int64,readdlm( MatchSubF[chk], skipstart=3 ) )
        am=convert.(Int64,readdlm( MatchMainF[chk], skipstart=3 ) )
        as = as[sp1,:]
        am = am[sp1,:]
    end

    #        @showprogress 2 "Loop over lensed catalog IDs..." for i=1:ng
    @time Threads.@threads for i=1:ng
        k = findall( id0t .== IDcat[i] )
        (length(k)<1) && continue

        ## fix ellipticities
        es1_cat[i] = e1t[ k[1] ]
        es2_cat[i] = e2t[ k[1] ]
        es =  es1_cat[i] + im*es2_cat[i]
        g = ( g1cat[i]  + im*g2cat[i] ) / (1.0-kappacat[i])
        if abs(g)<1.0
            e = ( es + g ) / ( 1.0 + conj(g)*es )
        else
            e = ( 1.0 + g*conj(es) ) / conj( es + g )
        end
        e1_cat[i] = real(e)
        e2_cat[i] = imag(e)
        if chk<4
            ## Match halos ellipticities
            kh = am[ k[1], 2 ] 
            ( (kh>0) && (kh<=nh) ) && ( Mh_main[ i ] = ms0[kh] )
            kh = as[ k[1], 2 ] 
            ( (kh>0) && (kh<=nh) ) && ( Mh_sub[ i ] = ms0[kh] )
        end
    end

    if born
        outf = joinpath( tmpd, root*"born_rescat_"*ch*".fits" )
    else
        outf = joinpath( tmpd, root*"rescat_"*ch*".fits" )
    end

    f = FITS( outf, "w" )
    data = Dict( "ID"=>IDcat, "RA"=>xcat, "DEC"=> ycat, "Z"=>zcat ,
                 "mstar"=>mscat, "Mh_sub"=>Mh_sub, "Mh_main"=>Mh_main,
                 "magnif"=>mucat, "kappa"=>kappacat,"gamma1"=>g1cat, "gamma2"=>g2cat,
                 "es1"=> es1_cat, "es2"=> es2_cat, "e1"=> e1_cat, "e2"=> e2_cat
                 )
    write( f, data )
    close( f )

end


#############################################
function read_res_cat( fn::AbstractString )
    f = FITS( fn )
    clis = ["ID","RA","DEC","Z","mstar","Mh_sub","Mh_main",
          "magnif","kappa","gamma1","gamma2","es1","es2","e1","e2" ]
    ncol = length(clis)
    aa = Vector{Any}(undef,ncol)
    for j = 1:ncol
        aa[j] = read(f[2],clis[j])
    end
    close(f)
    return aa,clis
end



#############################################
function write_csv_cat( ; of="summary_all.csv" )
    ff = open( joinpath( tmpd, of ), "w"  )
    clis=nothing
    for chk=1:length(chunck_list)
        ch = chunck_list[chk]
        (aa,clis) = read_res_cat( joinpath(tmpd,"rescat_"*ch*".fits") )
        ncol=length(clis)
        nr=length(aa[1])
        for i=1:nr
            s=@sprintf("%d,", (aa[1])[i] )
            for j=2:ncol-1
                s=s*@sprintf("%f,", (aa[j])[i])
            end
            s=s*@sprintf("%f", (aa[ncol])[i])
            @printf(ff,"%s\n",s)
        end
    end
    close(ff)
    return clis
end

#####################################################
function get_bbox2( x, y, W, cr1, cr2, cd1, cd2, nn; margin=0.02 )
    res = Vector{Int64}(undef,6)
    i0 = floor( Int64, Float32(x)/cd1 + cr1 )
    j0 = floor( Int64, Float32(y)/cd2 + cr2 )
    m = floor( Int64, margin/abs(cd1) )
    res[1] = max(m,i0-W)
    res[2] = min(nn[1]-m,i0+W)
    res[5] = max(m,min(i0,nn[1]-m))
    m = floor( Int64, margin/abs(cd2) )
    res[3] = max(m,j0-W)
    res[4] = min(nn[2]-m,j0+W)
    res[6] = max(m,min(j0,nn[2]-m))
    return res
end



##############################################
function rebin( arr::Array{Float32,2} ; factor::Int64=5 )
    n = size(arr)
    nn = (div(n[1],factor),div(n[2],factor))
    (nn1::Int64,nn2::Int64) = (nn[1],nn[2])
    res = Array{Float32,2}( undef, nn[1],nn[2] )
    ifsq::Float32 =1.0f0/factor^2
    @inbounds for j=1:nn2
        j1 = j*factor
        j0 = j1-factor+1
        @inbounds @simd for i=1:nn1
            i1 = i*factor
            i0 = i1-factor+1
            res[i,j] = ifsq * sum( view(arr,i0:i1,j0:j1) )
        end
    end
    return res
end

##############################################
function winc( i::Int64, j::Int64, flexf1::Array{Float32,2}, flexf2::Array{Float32,2}, flexg1::Array{Float32,2}, flexg2::Array{Float32,2}, g1::Array{Float32,2}, g2::Array{Float32,2}, kappa::Array{Float32,2},
               bb, cd2, acd1, rmin2, lrmin, dlr, nr, w::Vector{Float32}, ir::Vector{Int64}, kp::Vector{Bool} )
    w[1] = (j-bb[6])*cd2 #y
    w[2] = w[1]^2  #y2
    w[3] = (i-bb[5])*acd1 #x
    w[4] = w[3]^2  #x2
    w[5] = w[4] + w[2]  #r2
    r = sqrt(w[5]) 
    r3 = r^3
    x3 = w[3]^3
    y3 = w[1]^3
    kp[1] = w[5]>rmin2
    kp[1] || return
    w[6] = ( -(w[4]-w[2])*g1[i,j] - 2.f0*w[1]*w[3]*g2[i,j] ) / w[5]
    w[7] = kappa[i,j]
    w[8] = (-flexf1[i,j]*w[3]-w[1]*flexf2[i,j])/r  #-Re(F)
    #w[8] = (flexf1[i,j]*w[1]-w[3]*flexf2[i,j])/r   #-Im(F)
    w[9] = -flexg1[i,j]*(4*(x3/r3)-3*w[3]/r)-flexg2[i,j]*(-4*(y3/r3)+3*w[1]/r)  #-Re(G)
    #w[9] = -flexg2[i,j]*(4*(x3/r3)-3*w[3]/r)+flexg1[i,j]*(-4*(y3/r3)+3*w[1]/r)  #-Im(G)
    ir[1] = ceil( Int64, (0.5f0*log(w[5])-lrmin)/dlr )
    kp[1] = (ir[1]>=1) & (ir[1]<=nr)
end


##############################################
#prendre en compte ouverture du lightcone dx=1.0! de base rmin = 0.02 rmax = 50
function comp_gs_corr( ra, dec, fn::AbstractString; rf=4, dx=1.0f0, rmin=0.05f0, rmax=5f0, nr=50 )

    ng=length(ra)
    @assert length(dec)==ng

    ######### Deflection and derived information
    @time (alpha,zz) = read_bin( fn )
    nn=size(alpha)
    @info "done reading deflection map!"
    #@time jac=alpha2jac( alpha, scale=dx )
    #@info "done computing 1st order derivatives!"
    #@time jac=jac2der(jac)
    @time jac = alpha2jac(alpha,scale=dx)
    alpha=nothing ; GC.gc(true)
    @info "done computing jacobian matrix map"
    @time begin
        kappa = jac2kappa(jac)
        gamma1 = jac2gamma1(jac)
        gamma2 = jac2gamma2(jac)
	F1 = jac2F1(jac)
	F2 = jac2F2(jac)
	G1 = jac2G1(jac)
	G2 = jac2G2(jac)
	
    end
    jac=nothing ; GC.gc(true)
    @info "done computing raw convergence, shear and flexion maps"
    nn=(nn[1]/rf, nn[2]/rf)
    @time begin 
        kappa=rebin(kappa,factor=rf)
        gamma1=rebin(gamma1,factor=rf)
        gamma2=rebin(gamma2,factor=rf)
        F1=rebin(F1,factor=rf)
        F2=rebin(F2,factor=rf)
        G1=rebin(G1,factor=rf)
        G2=rebin(G2,factor=rf)	
    end
    @info "done computing rebinned convergence & shear maps"
    cd1 = -dx/nn[1] ; acd1=abs(cd1)
    cr1 = (nn[1]+1)/2f0
    cd2 = dx/nn[2]
    cr2 = (nn[2]+1)/2f0
    
    ######### Rad bins
    lrmin=log(rmin/60.f0) ; lrmax=log(rmax/60.f0) ;rmin2=(rmin/60.f0)^2
    dlr=(lrmax-lrmin)/(nr)
    W=ceil(Int64,rmax/abs(cd2)/60.f0)
    rbin = exp.( range( lrmin, stop=lrmax,length=nr+1 ) ) ;

    Sres=zeros(5,nr,ng)

    #######
    Threads.@threads for ig in eachindex(ra)
        bb=get_bbox2( ra[ig], dec[ig], W, cr1, cr2, cd1, cd2, nn )
	#afficher bb 1,2,3,4
        w=zeros(Float32,9)
        lres=zeros(Float64,5,nr)
        ir=[1]
        kp=[false]

        @inbounds for j=bb[3]:bb[4]
            @inbounds for i=bb[1]:bb[2]
                winc( i, j, F1, F2, G1, G2, gamma1, gamma2, kappa, bb, cd2, acd1, rmin2, lrmin, dlr, nr, w, ir, kp )
                if kp[1]
		    lres[1,ir[1]] += Float64(w[9]) #storing G+
		    lres[2,ir[1]] += Float64(w[8]) #storing F+
                    lres[3,ir[1]] += Float64(w[7]) #storing kappa
                    lres[4,ir[1]] += Float64(w[6]) #storing gamma+
                    lres[5,ir[1]] += 1.0           #counting pixels
                end
            end
        end
        Sres[:,:,ig] .= lres
    end
    #    Sres = pmap( ig->loc_calc(ra[ig],dec[ig]) , 1:ng )
    #    res = sdata(Sres)
    return Sres, rbin, kappa
end
