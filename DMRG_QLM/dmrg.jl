
using LinearAlgebra, TensorOperations, KrylovKit
using JLD


include("qlm_mpo.jl")
include("func_dmrg.jl")
include("chess_operator.jl")
s=16
D= 20
N = 9

A = randmps(N, s, D);
coupling_interaction=-10.0

#M =  mpoqlm(N ; coupling=1.  );
#M =  mpoqlm_fixed_boundary_contitions(N ; coupling=1.  );
M =  mpoqlm_fixed_with_interaction(N ; coupling=coupling_interaction  );


E, A, F = dmrgconvergence!(A, M ; verbose = true);

N=N-1
println("$N  $D    $E  $coupling_interaction ")
N=N+1

sp = [0. 1.; 0. 0.]
sm = [0. 0.; 1. 0.]
sz = [1. 0.; 0. -1.]
u = [1. 0.; 0. 1.]
pp = [1. 0.; 0. 0.]
pm = [0. 0.; 0. 1.]

O=kron(kron(kron(u,u),sz),sz)

winding_number=measure1siteoperator(A,O)
winding_number= deleteat!(winding_number , N)
#println(winding_number)



println(sum(winding_number) )


chess_down=chess_operator_down(N)

chess_up=chess_operator_up(N)

down = measure_mpo!(A,chess_down)

up = measure_mpo!(A,chess_up)


<<<<<<< HEAD
#println("Number_of_Plaquettes $N coupling_interaction $coupling_interaction Bond_dimention $D  Energy_GS  $E winding_number $(sum(winding_number)) chess_down $down chess_up $up ")
println("$N                   $coupling_interaction                 $D                 $E     $(sum(winding_number))      $down     $up ")
=======
<<<<<<< HEAD
println(v-winding_number[N])
"""
=======
println(down)
>>>>>>> a39e62167dcb6854e54bd822078f4eed97e9e945
>>>>>>> 19f7d5b0b51329bbad2d9c707875cd8bb8f2c411
