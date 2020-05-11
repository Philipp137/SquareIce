
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
#M =  mpoqlm_with_interaction_chemical_potential_and_magneti_field(N ; coupling=coupling_interaction , mu=chemical_potential, theta=theta);



#E, A, F = dmrgconvergence!(A, M ; verbose = true);


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


EA_mpo=operator_sublattice_energy(N,"A")

EB_mpo=operator_sublattice_energy(N,"B")

Oflip_mpo = operator_Oflip(N)

# Measure
EA    = measure_mpo!(A,EA_mpo)
EB    = measure_mpo!(A,EB_mpo)
Oflip = measure_mpo!(A,Oflip_mpo)
#println(E)

println("Number_Plaquettes  coupling       chemical        theta           Bond_dimention       Energy_GS                        winding_number                 EA              EB              Oflip ")
println("$(N-1)                  $coupling_interaction            $chemical_potential             $theta             $D                   $E            $(real(sum(winding_number)))              $(real(EA))     $(real(EB))     $(real(Oflip))")
