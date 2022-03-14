
using LinearAlgebra, TensorOperations, KrylovKit

include("src/mpos.jl")
include("src/func_dmrg.jl")
include("src/observables.jl")
include("src/io.jl")

using Main.mpos: mpoqlm

function measure(;config_name="",Lx= 0 , BondD=0 , Ly=0,theta = 0.0, coupling=0.0, mu=0.0, prefix = "SQ")

     if config_name != ""
         dict  =load_configuration( config_name)
     else
         mu_vec = [0.0, parse(Float64, mu)]
         coupling = parse(Float64, coupling)
         dict = load_configuration( Lx, BondD, coupling; theta , mu=mu_vec, Ly=Ly, prefix=prefix)
     end

     A =dict["MPS"]
     Lx =dict["Lx"]
     D =dict["D"]
     lambda =dict["lambda"]
     theta  =dict["theta"]
     mu_vec =dict["mu"]
     Ly =dict["Ly"]
     if Ly == 2
         s = 16
     else
         println("Warning!!! Ly >2 not implemented")
     end
     mpo = mpoqlm(Lx; Ly=Ly, coupling=lambda , mu=mu_vec , theta= theta);
     measure_op!(lambda , mu_vec,  theta  , s, D, 10*D, A, mpo)
end



#measure(config_name=ARGS[1])
measure(Lx= ARGS[1] , BondD=ARGS[2] , Ly=ARGS[3], coupling=ARGS[4], mu=ARGS[5], prefix = ARGS[6])
