using LinearAlgebra, TensorOperations, KrylovKit



function operator_flipp(N::Int;   coupling=1. )
    sp = [0. 1.; 0. 0.]
    sm = [0. 0.; 1. 0.]
    sz = [1. 0.; 0. -1.]
    u = [1. 0.; 0. 1.]
    pp = [1. 0.; 0. 0.]
    pm = [0. 0.; 0. 1.]
    
    s=16
    D = 6
    M = zeros(D, s, D, s) ;
    MN = zeros(D, s, D, s) ;

    I= operator_2=kron(kron(kron(u,u),u),u)


    M[1,:,1,:] = I ; M[ D,:,D,:] = I
    MN[1,:,1,:] = I ;    MN[ D,:,D,:] = I


    operator_5=kron(kron(kron(pm,u),pm),pp)
    operator_6=kron(kron(kron(pp,u),u),u)
    
    operator_7=-1.0 * kron(kron(kron(u,pm),pp),pm)
    operator_8=kron(kron(kron(u,pp),u),u)
    

    operator_5_dag=kron(kron(kron(pp,u),pp),pm)
    operator_6_dag=kron(kron(kron(pm,u),u),u)
    
    operator_7_dag=-1.0*kron(kron(kron(u,pp),pm),pp)
    operator_8_dag=kron(kron(kron(u,pm),u),u)

    M[1,:,2,:] = operator_5
    M[2,:,D,:] = operator_6

    M[1,:,3,:] = operator_7
    M[3,:,D,:] = operator_8

    M[1,:,4,:] = operator_5_dag
    M[4,:,D,:] = operator_6_dag

    M[1,:,5,:] = operator_7_dag
    M[5,:,D,:] = operator_8_dag


    MN[2,:,D,:]  = operator_6
    MN[3,:,D,:] = operator_8
    MN[4,:,D,:] = operator_6_dag
    MN[5,:,D,:] = operator_8_dag

    mpo= [ M[1:1,:,:,:] ]

    for site = 2:N-1;
        push!(mpo, (-1.0)^(site)*M[:,:,:,:])
    end

    push!(mpo, (-1.0)^(N)*MN[:,:,D:D,:])

    return mpo
end
