
using LinearAlgebra, TensorOperations, KrylovKit

include("func_dmrg.jl")
include("mpo_gn.jl")


function measure_chiral_condensate(A)
    N = length(A)
    sz = [1.0 0.0; 0.0 -1.0]

    expval = measure1siteoperator(A, sz)

    chiral_condensate=0.0
    
    mult=-1.0
    
    for b = 1:N;
        chiral_condensate+= mult*expval[b]
        mult=-1*mult
    end

    return chiral_condensate
end



function measure_operator_mpo(A, H)
    N = length(A)

    expval = 0.

    F = Vector{Any}(undef, N+2)
    F[1] = fill!(similar(H[1], (1,1,1)), 1)
    F[N+2] = fill!(similar(H[1], (1,1,1)), 1)
    for k = N:-1:1
        F[k+1] = updaterightenv(A[k], H[k], F[k+2])
    end
    
    @tensor v = F[1][α,a,α']*F[2][α,a,α']
  
    return v
end

function measure_interaction_term(A,O)
    N = length(A)
    interaction=measure_operator_mpo(A, O)
    return interaction
end


function measure_id(A)
    N = length(A)
    sz = [1.0 0.0; 0.0 1.0]

    expval = measure1siteoperator(A, sz)

    chiral_condensate=0.0
    
   
    for b = 1:N;
        chiral_condensate+= expval[b]
    end


    return chiral_condensate/N
end