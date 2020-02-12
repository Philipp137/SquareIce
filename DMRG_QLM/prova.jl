
using LinearAlgebra, TensorOperations, KrylovKit



include("qlm_mpo.jl")
include("func_dmrg.jl")
include("chess_operator.jl")
s=16
D= 20
N = 13
D_max=100

#A = randmps(N, s, D);
A = randmps(N, s, D, Complex{Float64});
coupling_interaction=-1.0
chemical_potential=0.0
theta=0.0

#M =  mpoqlm(N ; coupling=1.  );
#M =  mpoqlm_fixed_boundary_contitions(N ; coupling=1.  );
#M =  mpoqlm_fixed_with_interaction(N ; coupling=coupling_interaction  );
#M =  mpoqlm_with_interaction_and_chemical_potential(N ; coupling=coupling_interaction  , mu=chemical_potential=1.);
#M =  mpoqlm_with_interaction_and_chemical_potential(N ; coupling=coupling_interaction , mu=chemical_potential);
M =  mpoqlm_with_interaction_chemical_potential_and_magneti_field(N ; coupling=coupling_interaction , mu=chemical_potential, theta=theta);

dmrgconvergence_in_D!(D, D_max , A, M )
