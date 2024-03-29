
using LinearAlgebra, TensorOperations, KrylovKit



include("qlm_mpo.jl")
include("func_dmrg.jl")
include("chess_operator.jl")
s=16
D= 20
N = 5

#A = randmps(N, s, D);
A = randmps(N, s, D, Complex{Float64});
coupling_interaction=1.0
chemical_potential=0.0
theta=0.0

#M =  mpoqlm(N ; coupling=1.  );
#M =  mpoqlm_fixed_boundary_contitions(N ; coupling=1.  );
#M =  mpoqlm_fixed_with_interaction(N ; coupling=coupling_interaction  );
#M =  mpoqlm_with_interaction_and_chemical_potential(N ; coupling=coupling_interaction  , mu=chemical_potential=1.);
#M =  mpoqlm_with_interaction_and_chemical_potential(N ; coupling=coupling_interaction , mu=chemical_potential);
M =  mpoqlm_with_interaction_chemical_potential_and_magneti_field(N ; coupling=coupling_interaction , mu=chemical_potential, theta=theta);



E, A, F = dmrgconvergence!(A, M ; verbose = true);


N=N-1
#println("$N  $D    $E  $coupling_interaction ")
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



#println(sum(winding_number) )


chess_down=chess_operator_down(N)

chess_up=chess_operator_up(N)

down = measure_mpo!(A,chess_down)

up = measure_mpo!(A,chess_up)
#println(E)

println("Number_Plaquettes  coupling       chemical        theta           Bond_dimention       Energy_GS                        winding_number                 chess_down              chess_up ")
println("$(N-1)                  $coupling_interaction            $chemical_potential             $theta             $D                   $E            $(real(sum(winding_number)))              $(real(down))     $(real(up)) ")
