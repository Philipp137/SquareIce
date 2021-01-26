
using LinearAlgebra, TensorOperations, KrylovKit

include("src/mpos.jl")
include("src/func_dmrg.jl")
include("src/observables.jl")
include("src/io.jl")

using Main.mpos: mpoqlm

function dostuff(N ,  D , coupling; storage="./data/")
    s = 16
    D_max = 10*D
    A = randmps(N, s, D);
    chemical_potential = 0.0
    theta = 0.0
    mpo = mpoqlm(N, coupling=coupling , mu=chemical_potential , theta= theta);
    #mpo = mpoqlm(N, coupling=coupling );
    E, A, F, counter =   dmrgconvergence_in_D!(s, D, D_max, A, mpo);
    E, A, F , counter = dmrgconvergence_in_D_and_measure_op!(coupling, chemical_potential,  theta, s, D, D_max, A, mpo);
    save_configuration(A, N , D, coupling, theta= theta, mu=chemical_potential, prefix =string(storage,"/SQ"))

end



dostuff(parse(Int, ARGS[1]) , parse(Int, ARGS[2]) , parse(Float64, ARGS[3]))
