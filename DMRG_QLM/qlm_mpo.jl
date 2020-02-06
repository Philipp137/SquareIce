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

function mpoqlm_fixed_with_interaction(N::Int;   coupling=-1. )
    sp = [0. 1.; 0. 0.]
    sm = [0. 0.; 1. 0.]
    sz = [1. 0.; 0. -1.]
    u = [1. 0.; 0. 1.]
    pp = [1. 0.; 0. 0.]
    pm = [0. 0.; 0. 1.]
    
    s=16
    D = 10
    M = zeros(D, s, D, s) ;
    MN = zeros(D, s, D, s) ;

    I= operator_2=kron(kron(kron(u,u),u),u)


    M[1,:,1,:] = I ;    M[ D,:,D,:] = I
    MN[1,:,1,:] = I ;    MN[ D,:,D,:] = I

    operator_1=-1*kron(kron(kron(sm,u),sm),sp)
    operator_2=kron(kron(kron(sp,u),u),u)
    operator_3=-1*kron(kron(kron(u,sm),sp),sm)
    operator_4=kron(kron(kron(u,sp),u),u)

    operator_5=coupling*kron(kron(kron(pm,u),pm),pp)
    operator_6=kron(kron(kron(pp,u),u),u)
    
    operator_7=coupling*kron(kron(kron(u,pm),pp),pm)
    operator_8=kron(kron(kron(u,pp),u),u)
    

    operator_1_dag=-1*kron(kron(kron(sp,u),sp),sm)    
    operator_2_dag=kron(kron(kron(sm,u),u),u)
    
    operator_3_dag=-1*kron(kron(kron(u,sp),sm),sp)
    operator_4_dag=kron(kron(kron(u,sm),u),u)

    operator_5_dag=coupling*kron(kron(kron(pp,u),pp),pm)
    operator_6_dag=kron(kron(kron(pm,u),u),u)
    
    operator_7_dag=coupling*kron(kron(kron(u,pp),pm),pp)
    operator_8_dag=kron(kron(kron(u,pm),u),u)

    projector_N=kron(kron(kron(u,u),pp),pm)


    M[1,:,2,:] = operator_1
    M[2,:,D,:] = operator_2

    M[1,:,3,:] = operator_3
    M[3,:,D,:] = operator_4

    M[1,:,4,:] = operator_1_dag
    M[4,:,D,:] = operator_2_dag

    M[1,:,5,:] = operator_3_dag
    M[5,:,D,:] = operator_4_dag

    M[1,:,6,:] = operator_5
    M[6,:,D,:] = operator_6

    M[1,:,7,:] = operator_7
    M[7,:,D,:] = operator_8

    M[1,:,8,:] = operator_5_dag
    M[8,:,D,:] = operator_6_dag

    M[1,:,9,:] = operator_7_dag
    M[9,:,D,:] = operator_8_dag

    MN[2,:,D,:] = operator_2
    MN[3,:,D,:] = operator_4
    MN[1,:,D,:] = projector_N
    MN[4,:,D,:] = operator_2_dag
    MN[5,:,D,:] = operator_4_dag
    MN[6,:,D,:] = operator_6
    MN[7,:,D,:] = operator_8
    MN[8,:,D,:] = operator_6_dag
    MN[9,:,D,:] = operator_8_dag

    mpo= [ M[1:1,:,:,:] ]

    for site = 2:N-1;
        push!(mpo, M[:,:,:,:])
    end

    push!(mpo, MN[:,:,D:D,:])

    return mpo
end

function mpoqlm_with_interaction_and_chemical_potential(N::Int;   coupling=-1. , mu=-1 )
    sp = [0. 1.; 0. 0.]
    sm = [0. 0.; 1. 0.]
    sz = [1. 0.; 0. -1.]
    u = [1. 0.; 0. 1.]
    pp = [1. 0.; 0. 0.]
    pm = [0. 0.; 0. 1.]
    
    s=16
    D = 10
    M = zeros(D, s, D, s) ;
    MN = zeros(D, s, D, s) ;

    I= operator_2=kron(kron(kron(u,u),u),u)
    mu_1=kron(kron(kron(sz,u),u),u)
    mu_2=kron(kron(kron(u,sz),u),u)
    mu_3=kron(kron(kron(u,u),sz),u)
    mu_4=kron(kron(kron(u,u),u),sz)

    M[1,:,1,:] = I ;    M[ D,:,D,:] = I
    MN[1,:,1,:] = I ;    MN[ D,:,D,:] = I

    operator_1=-1*kron(kron(kron(sm,u),sm),sp)
    operator_2=kron(kron(kron(sp,u),u),u)
    operator_3=-1*kron(kron(kron(u,sm),sp),sm)
    operator_4=kron(kron(kron(u,sp),u),u)

    operator_5=coupling*kron(kron(kron(pm,u),pm),pp)
    operator_6=kron(kron(kron(pp,u),u),u)
    
    operator_7=coupling*kron(kron(kron(u,pm),pp),pm)
    operator_8=kron(kron(kron(u,pp),u),u)
    

    operator_1_dag=-1*kron(kron(kron(sp,u),sp),sm)    
    operator_2_dag=kron(kron(kron(sm,u),u),u)
    
    operator_3_dag=-1*kron(kron(kron(u,sp),sm),sp)
    operator_4_dag=kron(kron(kron(u,sm),u),u)

    operator_5_dag=coupling*kron(kron(kron(pp,u),pp),pm)
    operator_6_dag=kron(kron(kron(pm,u),u),u)
    
    operator_7_dag=coupling*kron(kron(kron(u,pp),pm),pp)
    operator_8_dag=kron(kron(kron(u,pm),u),u)

    projector_N=kron(kron(kron(u,u),pp),pm)


    M[1,:,D,:] = mu*(mu_1+mu_2+mu_3+mu_4)

    M[1,:,2,:] = operator_1
    M[2,:,D,:] = operator_2

    M[1,:,3,:] = operator_3
    M[3,:,D,:] = operator_4

    M[1,:,4,:] = operator_1_dag
    M[4,:,D,:] = operator_2_dag

    M[1,:,5,:] = operator_3_dag
    M[5,:,D,:] = operator_4_dag

    M[1,:,6,:] = operator_5
    M[6,:,D,:] = operator_6

    M[1,:,7,:] = operator_7
    M[7,:,D,:] = operator_8

    M[1,:,8,:] = operator_5_dag
    M[8,:,D,:] = operator_6_dag

    M[1,:,9,:] = operator_7_dag
    M[9,:,D,:] = operator_8_dag

    MN[1,:,D,:] = mu*(mu_1+mu_2)

    MN[2,:,D,:] = operator_2
    MN[3,:,D,:] = operator_4
    MN[1,:,D,:] = projector_N
    MN[4,:,D,:] = operator_2_dag
    MN[5,:,D,:] = operator_4_dag
    MN[6,:,D,:] = operator_6
    MN[7,:,D,:] = operator_8
    MN[8,:,D,:] = operator_6_dag
    MN[9,:,D,:] = operator_8_dag

    mpo= [ M[1:1,:,:,:] ]

    for site = 2:N-1;
        push!(mpo, M[:,:,:,:])
    end

    push!(mpo, MN[:,:,D:D,:])

    return mpo
end
