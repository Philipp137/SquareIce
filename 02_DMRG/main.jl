
using LinearAlgebra, TensorOperations, KrylovKit

include("src/mpos.jl")
include("src/func_dmrg.jl")
include("src/observables.jl")
include("src/io.jl")

using Main.mpos: mpoqlm

function dostuff(Lx , D , Ly, coupling, mu; storage="./data/")

    Nlinks = Ly * 2
    Nstates = 2^(Nlinks) # number of states for one chain element

    D_max = 10*D
    A = randmps(Lx, Nstates, D);
    chemical_potential = [0, mu]
    theta = 0.0
    mpo = mpoqlm( Lx, Ly = Ly, coupling=coupling , mu=chemical_potential , theta= theta);
    #mpo = mpoqlm(N, coupling=coupling );
    E, A, F, counter =   dmrgconvergence_in_D!(Nstates, D, D_max,  A, mpo, Ly = Ly);
    E, A, F , counter = dmrgconvergence_in_D_and_measure_op!(coupling, chemical_potential,  theta, Nstates, D, D_max, A, mpo, Ly = Ly);
    save_configuration(A, Lx , D, coupling, theta= theta, mu=chemical_potential, Ly = Ly, prefix =string(storage,"/SQ"))

end



dostuff(parse(Int, ARGS[1]) , parse(Int, ARGS[2]) , parse(Int, ARGS[3]) , parse(Float64, ARGS[4]), parse(Float64,ARGS[5]))
