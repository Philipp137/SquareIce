
using LinearAlgebra, TensorOperations, KrylovKit

include("qlm_mpo.jl")
include("func_dmrg.jl")
include("chess_operator.jl")
function dostuff()

    s = 16
    D = 40
    N = 101
    D_max = 101
    coupling_interaction = -5.0

    println("Number_Plaquettes  coupling       chemical        theta           Bond_dimention       Energy_GS                        winding_number                 chess_down              chess_up ")

    while  coupling_interaction < -4.0
        coupling = coupling_interaction
        A = randmps(N, s, D);
        chemical_potential = 0.0
        theta = 0.0
        M =  mpoqlm_fixed_with_interaction(N ; coupling = coupling);
        E, A, F, counter =   dmrgconvergence_in_D!(s, D, D_max, A, M)
        dmrgconvergence_in_D_and_measure_op!(coupling_interaction, chemical_potential,  theta, s, D, D_max, A, M)
        coupling_interaction = coupling_interaction + 7.0
            
    end

end

dostuff()
