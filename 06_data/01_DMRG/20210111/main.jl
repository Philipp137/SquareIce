using LinearAlgebra, TensorOperations, KrylovKit

include("../../../02_DMRG/mpos.jl")
#include("../../../02_DMRG/archiv/qlm_mpo.jl")
include("../../../02_DMRG/func_dmrg.jl")
#include("./scan_ft_SI/chess_operator.jl")
include("../../../02_DMRG/observables.jl")

using Main.mpos: mpoqlm

function dostuff(N ,  D , coupling)
    s = 16
    D_max = 10*D
    A = randmps(N, s, D);
    chemical_potential = 0.0
    theta = 0.0
    mpo = mpoqlm(N, coupling=coupling , mu=chemical_potential , theta= theta);
    #mpo = mpoqlm(N, coupling=coupling );
    E, A, F, counter =   dmrgconvergence_in_D!(s, D, D_max, A, mpo);
    dmrgconvergence_in_D_and_measure_op!(coupling, chemical_potential,  theta, s, D, D_max, A, mpo);
end



dostuff(parse(Int, ARGS[1]) , parse(Int, ARGS[2]) , parse(Float64, ARGS[3]))
