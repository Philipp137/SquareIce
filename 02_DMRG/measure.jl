
using LinearAlgebra, TensorOperations, KrylovKit

include("src/mpos.jl")
include("src/func_dmrg.jl")
include("src/observables.jl")
include("src/io.jl")

using Main.mpos: mpoqlm

function measure(config_name="./data/my_config.jld2")

     dict  =load_configuration( config_name)
     A =dict["MPS"]
     Lx =dict["Lx"]
     D =dict["D"]
     lambda =dict["lambda"]
     theta  =dict["theta"]
     mu =dict["mu"]
     Ly =dict["Ly"]
     if Ly == 2
         s = 16
     else
         @printf("Warning!!! Ly >2 not implemented")
     end
     mpo = mpoqlm(Lx, coupling=lambda , mu=mu , theta= theta);
     measure_op!(lambda ,mu  ,  theta  , s, D, 10*D, A, mpo)
end



measure(ARGS[1])
