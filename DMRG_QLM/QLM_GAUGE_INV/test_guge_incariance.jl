using LinearAlgebra
using Arpack
using SparseArrays
using KrylovKit

function projector()
    sz = [1. 0.; 0. -1.]
    u = [1. 0.; 0. 1.]
#    projector_low= kron(kron(kron(sz,-sz),sz),-sz) 
    projector_low= 0.5*( kron(kron(kron(kron(kron(kron( kron(u,u),-sz),u),u),u),u),u) +  kron(kron(kron(kron(kron(kron( kron(u,u),u),u),sz),u),u),u) +   kron(kron(kron(kron(kron(kron( kron(u,u),u),u),u),u),sz),u)  + kron(kron(kron(kron(kron(kron( kron(u,u),u),u),u),-sz),u),u) )


    
    return projector_low
end


A=eigvals(projector())
B=eigvecs(projector())

len=length(A)

C = Float64[]

for i  in 1:len
    if abs( A[i])< 0.1 

#        println(i)
        push!(C,B[i] )

    end

end

println(C  )