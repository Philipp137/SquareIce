
using LinearAlgebra, TensorOperations, KrylovKit


function mpoqlm(N::Int;   coupling=1. )
    sp = [0. 1.; 0. 0.]
    sm = [0. 0.; 1. 0.]
    sz = [1. 0.; 0. -1.]
    u = [1. 0.; 0. 1.]

    s=16
    D = 6
    M = zeros(D, s, D, s) ;

    I= operator_2=kron(kron(kron(u,u),u),u)


    M[1,:,1,:] = I ;    M[ D,:,D,:] = I

    operator_1=kron(kron(kron(sm,u),sm),sp)
    operator_2=kron(kron(kron(sp,u),u),u)
    operator_3=kron(kron(kron(u,sm),sp),sm)
    operator_4=kron(kron(kron(u,sp),u),u)
    operator_1_dag=kron(kron(kron(sp,u),sp),sm)    
    operator_2_dag=kron(kron(kron(sm,u),u),u)
    operator_3_dag=kron(kron(kron(u,sp),sm),sp)
    operator_4_dag=kron(kron(kron(u,sm),u),u)

    M[1,:,2,:] = operator_1
    M[2,:,D,:] = operator_2

    M[1,:,3,:] = operator_3
    M[3,:,D,:] = operator_4

    M[1,:,4,:] = operator_1_dag
    M[4,:,D,:] = operator_2_dag

    M[1,:,5,:] = operator_3_dag
    M[5,:,D,:] = operator_4_dag

    mpo= [ M[1:1,:,:,:] ]

    for site = 2:N-1;
        push!(mpo, M[:,:,:,:])
    end

    push!(mpo, M[:,:,D:D,:])

    return mpo
end

function mpoqlm_fixed_boundary_contitions(N::Int;   coupling=1. )
    sp = [0. 1.; 0. 0.]
    sm = [0. 0.; 1. 0.]
    sz = [1. 0.; 0. -1.]
    u = [1. 0.; 0. 1.]
    pp = [1. 0.; 0. 0.]
    pm = [0. 0.; 0. 1.]
    
    s=16
    D = 6
    M = zeros(D, s, D, s) ;
    M1 = zeros(D, s, D, s) ;
    MN = zeros(D, s, D, s) ;

    I= operator_2=kron(kron(kron(u,u),u),u)


    M[1,:,1,:] = I ;    M[ D,:,D,:] = I
    M1[1,:,1,:] = I ;    M1[ D,:,D,:] = I
    MN[1,:,1,:] = I ;    MN[ D,:,D,:] = I

    operator_1=kron(kron(kron(sm,u),sm),sp)
    operator_2=kron(kron(kron(sp,u),u),u)
    operator_3=kron(kron(kron(u,sm),sp),sm)
    operator_4=kron(kron(kron(u,sp),u),u)
    operator_1_dag=kron(kron(kron(sp,u),sp),sm)    
    operator_2_dag=kron(kron(kron(sm,u),u),u)
    operator_3_dag=kron(kron(kron(u,sp),sm),sp)
    operator_4_dag=kron(kron(kron(u,sm),u),u)

    projector_1=kron(kron(kron(u,u),pp),pm)
    projector_N=kron(kron(kron(u,u),pp),pm)


    M1[1,:,D,:] = projector_1

    M[1,:,2,:] = operator_1
    M[2,:,D,:] = operator_2

    M[1,:,3,:] = operator_3
    M[3,:,D,:] = operator_4

    M[1,:,4,:] = operator_1_dag
    M[4,:,D,:] = operator_2_dag

    M[1,:,5,:] = operator_3_dag
    M[5,:,D,:] = operator_4_dag

    MN[2,:,D,:] = operator_2
    MN[3,:,D,:] = operator_4
    MN[1,:,D,:] = projector_N
    MN[4,:,D,:] = operator_2_dag
    MN[5,:,D,:] = operator_4_dag

    mpo= [ M[1:1,:,:,:] ]

    for site = 2:N-1;
        push!(mpo, M[:,:,:,:])
    end

    push!(mpo, MN[:,:,D:D,:])

    return mpo
end
