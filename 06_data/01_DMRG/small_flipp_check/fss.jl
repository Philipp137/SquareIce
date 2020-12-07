using LinearAlgebra, TensorOperations, KrylovKit

include("./scan_ft_SI/qlm_mpo.jl")
include("./scan_ft_SI/func_dmrg.jl")
#include("./scan_ft_SI/chess_operator.jl")
include("./scan_ft_SI/observables.jl")
function dostuff(N ,  D , coupling)
    s = 16
    D_max = 10*D
    println("Number_Plaquettes  coupling       chemical        theta           Bond_dimention       Energy_GS                        winding_number                 flipp ")
    A = randmps(N, s, D);
    chemical_potential = 0.0
    theta = 0.0
    M =  mpoqlm_fixed_with_interaction(N ; coupling = coupling);
    E, A, F, counter =   dmrgconvergence_in_D!(s, D, D_max, A, M)
    dmrgconvergence_in_D_and_measure_op!(coupling, chemical_potential,  theta, s, D, D_max, A, M)

end



dostuff(parse(Int, ARGS[1]) , parse(Int, ARGS[2]) , parse(Float64, ARGS[3]))
