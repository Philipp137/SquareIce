
using LinearAlgebra, TensorOperations, KrylovKit
using JLD


include("qlm_mpo.jl")
include("func_dmrg.jl")

s=10
D= 200
N = 20

A = randmps(N, s, D);


M =  mpoqlm_gauge_inv(N ; coupling=1.  );
#M = mpogn_continuous_mu(N; mass=mass,  coupling_V=coupling_V,  coupling_A=coupling_A, mu=mu)


E, A, F = dmrgconvergence!(A, M ; verbose = true);


#dt=0.01
#A, F = tdvp1sweep!(dt, A, M; verbose = false);

println("$N  $D    $E   ")







