
using LinearAlgebra, TensorOperations, KrylovKit

include("src/mpos.jl")
include("src/func_dmrg.jl")
include("src/observables.jl")
include("src/io.jl")

using Main.mpos: mpoqlm

function measure(;config_name="./data/my_config.jld2",Lx= 0 , BondD=0 , Ly=0,theta = 0, coupling=0, mu=0, prefix = "SQ")

     if config_name != ""
         dict  =load_configuration( config_name)
     else
         dict = load_configuration( Lx, BondD, coupling; theta , mu, Ly, prefix)
     end

     A =dict["MPS"]
     Lx =dict["Lx"]
     D =dict["D"]
     lambda =dict["lambda"]
     theta  =dict["theta"]
     mu =dict["mu"]
     Ly =dict["Ly"]
     chemical_potential = [0, mu]
     if Ly == 2
         s = 16
     else
         println("Warning!!! Ly >2 not implemented")
     end
     mpo = mpoqlm(Lx; Ly=Ly, coupling=lambda , mu=chemical_potential , theta= theta);
     measure_op!(lambda , chemical_potential,  theta  , s, D, 10*D, A, mpo)
end



#measure(config_name=ARGS[1])
measure(Lx= ARGS[1] , BondD=ARGS[2] , Ly=ARGS[3], coupling=ARGS[4], mu=ARGS[5], prefix = ARGS[7])
