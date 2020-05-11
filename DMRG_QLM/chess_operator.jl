using LinearAlgebra, TensorOperations, KrylovKit



function chess_operator_down(N::Int )
    sp = [0. 1.; 0. 0.]
    sm = [0. 0.; 1. 0.]
    sz = [1. 0.; 0. -1.]
    u = [1. 0.; 0. 1.]
    pp = [1. 0.; 0. 0.]
    pm = [0. 0.; 0. 1.]

    s=16
    D = 6
    M = zeros(D, s, D, s) ;
    Meven = zeros(D, s, D, s) ;
    Modd = zeros(D, s, D, s) ;
    M1 = zeros(D, s, D, s) ;
    MN = zeros(D, s, D, s) ;

    I= operator_2=kron(kron(kron(u,u),u),u)


    Meven[1,:,1,:] = I ; Meven[ D,:,D,:] = I
    Modd[1,:,1,:] = I ;     Modd[ D,:,D,:] = I
    MN[1,:,1,:] = I ;    MN[ D,:,D,:] = I


    operator_5=kron(kron(kron(pm,u),pm),pp)
    operator_6=kron(kron(kron(pp,u),u),u)

    operator_7=kron(kron(kron(u,pm),pp),pm)
    operator_8=kron(kron(kron(u,pp),u),u)


    operator_5_dag=kron(kron(kron(pp,u),pp),pm)
    operator_6_dag=kron(kron(kron(pm,u),u),u)

    operator_7_dag=kron(kron(kron(u,pp),pm),pp)
    operator_8_dag=kron(kron(kron(u,pm),u),u)

    Meven[1,:,2,:] = operator_5
    Modd[2,:,D,:]  = operator_6

    Modd[1,:,3,:] = operator_7
    Meven[3,:,D,:] = operator_8

    Meven[1,:,4,:] = operator_5_dag
    Modd[4,:,D,:] =  -1.*operator_6_dag

    Modd[1,:,5,:] = operator_7_dag
    Meven[5,:,D,:] = -1.*operator_8_dag


    MN[2,:,D,:]  = operator_6

    MN[3,:,D,:] = operator_8

    MN[4,:,D,:] = -1.*operator_6_dag

    MN[5,:,D,:] = -1.*operator_8_dag
    mpo= [ Modd[1:1,:,:,:] ]

    for site = 2:N-1;
        if iseven(site)
            push!(mpo, -1*Meven[:,:,:,:])
        else
            push!(mpo, Modd[:,:,:,:])
        end
    end

    push!(mpo, MN[:,:,D:D,:])

    return mpo
end


function chess_operator_up(N::Int)
    sp = [0. 1.; 0. 0.]
    sm = [0. 0.; 1. 0.]
    sz = [1. 0.; 0. -1.]
    u = [1. 0.; 0. 1.]
    pp = [1. 0.; 0. 0.]
    pm = [0. 0.; 0. 1.]

    s=16
    D = 6
    M = zeros(D, s, D, s) ;
    Meven = zeros(D, s, D, s) ;
    Modd = zeros(D, s, D, s) ;
    M1 = zeros(D, s, D, s) ;
    MN = zeros(D, s, D, s) ;

    I= operator_2=kron(kron(kron(u,u),u),u)


    Meven[1,:,1,:] = I ; Meven[ D,:,D,:] = I
    Modd[1,:,1,:] = I ;     Modd[ D,:,D,:] = I
    MN[1,:,1,:] = I ;    MN[ D,:,D,:] = I


    operator_5=kron(kron(kron(pm,u),pm),pp)
    operator_6=kron(kron(kron(pp,u),u),u)

    operator_7=kron(kron(kron(u,pm),pp),pm)
    operator_8=kron(kron(kron(u,pp),u),u)


    operator_5_dag=kron(kron(kron(pp,u),pp),pm)
    operator_6_dag=kron(kron(kron(pm,u),u),u)

    operator_7_dag=kron(kron(kron(u,pp),pm),pp)
    operator_8_dag=kron(kron(kron(u,pm),u),u)

    Meven[1,:,2,:] = operator_5
    Modd[2,:,D,:]  = operator_6

    Modd[1,:,3,:] = operator_7
    Meven[3,:,D,:] = operator_8

    Meven[1,:,4,:] = operator_5_dag
    Modd[4,:,D,:] =  -1.*operator_6_dag

    Modd[1,:,5,:] = operator_7_dag
    Meven[5,:,D,:] = -1.*operator_8_dag


    MN[2,:,D,:]  = operator_6

    MN[3,:,D,:] = operator_8

    MN[4,:,D,:] = -1.*operator_6_dag

    MN[5,:,D,:] = -1.*operator_8_dag

    #even and odd are inveted in the computation of the mpo since it has the same method before

    mpo= [ Meven[1:1,:,:,:] ]

    for site = 2:N-1;
        if iseven(site)
            push!(mpo, Modd[:,:,:,:])
        else
            push!(mpo, Meven[:,:,:,:])
        end
    end

    push!(mpo, MN[:,:,D:D,:])

    return mpo
end
