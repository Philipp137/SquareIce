
using LinearAlgebra, TensorOperations, KrylovKit
using JLD


include("qlm_mpo.jl")
include("func_dmrg.jl")

s=16
D= 10
N = 7

A = randmps(N, s, D);


#M =  mpoqlm(N ; coupling=1.  );
#M =  mpoqlm_fixed_boundary_contitions(N ; coupling=1.  );
M =  mpoqlm_fixed_with_interaction(N ; coupling=-0.0  );


E, A, F = dmrgconvergence!(A, M ; verbose = true);

N=N-1
println("$N  $D    $E   ")
N=N+1

sp = [0. 1.; 0. 0.]
sm = [0. 0.; 1. 0.]
sz = [1. 0.; 0. -1.]
u = [1. 0.; 0. 1.]
pp = [1. 0.; 0. 0.]
pm = [0. 0.; 0. 1.]

O=kron(kron(kron(u,u),sz),sz)

winding_number=measure1siteoperator(A,O)

println(winding_number)



println(sum(winding_number) )

E = measure_mpo!(A,M)

println(E)