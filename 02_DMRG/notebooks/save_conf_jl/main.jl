using JLD
using LinearAlgebra, TensorOperations, KrylovKit



include("/Users/paolostornati/Phd/SquareIce/02_DMRG/func_dmrg.jl");




N= 10
s = 16
D = 100

A = randmps(N, s, D, Complex{Float64});

lambda= 0.5
theta = 0.4

save("./storage$N.jld", "N", N, "MPS", A , "D" , D , "lambda" , lambda , "theta" , theta , "s", s)

d = load("./storage$N.jld")

println( get(d, "MPS", true) )

