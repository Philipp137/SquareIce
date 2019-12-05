
using LinearAlgebra, TensorOperations, KrylovKit
using JLD


include("qlm_mpo.jl")
include("func_dmrg.jl")

s=16
D= 50
N = 9

A = randmps(N, s, D);


#M =  mpoqlm(N ; coupling=1.  );
M =  mpoqlm_fixed_boundary_contitions(N ; coupling=1.  );


E, A, F = dmrgconvergence!(A, M ; verbose = true);


println("$N  $D    $E   ")




