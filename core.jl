using Libdl
#const sllib=joinpath(homedir(),"lib","libsoftlens."*Libdl.dlext)
const sllib=joinpath(ENV["SOFTLENS_DIR"],"lib","libsoftlens."*Libdl.dlext)

const sl_version=(s="_"^40 ; ccall( (:jsl_version, sllib), Cvoid, (Ptr{UInt8}, Ref{Int32}),
           s, Int32(length(s)) ) ; strip(s) )
