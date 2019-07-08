
using LinearAlgebra, TensorOperations, KrylovKit
using JLD


include("qlm_mpo.jl")
include("func_dmrg.jl")

s=16
D= 20
N = 100

A = randmps(N, s, D);


M =  mpoqlm(N ; coupling=1.  );
#M = mpogn_continuous_mu(N; mass=mass,  coupling_V=coupling_V,  coupling_A=coupling_A, mu=mu)


E, A, F = dmrgconvergence!(A, M ; verbose = true);


println("$N  $D    $E   ")




