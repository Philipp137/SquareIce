using LinearAlgebra, TensorOperations, KrylovKit



function operator_flipp(N::Int;   coupling=1. )

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



# 1. Oflip = < \sum_plaq f_\plaq^2 >;
#    with f_\plaq = 1 if the plaquette is flippable, and 0 if not.
function operator_Oflip(N::Int )

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

    operator_7=kron(kron(kron(u,pm),pp),pm)
    operator_8=kron(kron(kron(u,pp),u),u)


    operator_5_dag=kron(kron(kron(pp,u),pp),pm)
    operator_6_dag=kron(kron(kron(pm,u),u),u)

    operator_7_dag=kron(kron(kron(u,pp),pm),pp)
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
        push!(mpo, M[:,:,:,:])
    end

    push!(mpo, MN[:,:,D:D,:])

    return mpo
end


# 1.eA and eB is the same as Oflip but measured only on a single sub-lattice.
function operator_sublattice_energy(N::Int, sublattice="A" )

    u = [1. 0.; 0. 1.]
    pp = [1. 0.; 0. 0.]
    pm = [0. 0.; 0. 1.]

    s=16
    D = 4
    M1 = zeros(D, s, D, s) ;
    MN1 = zeros(D, s, D, s) ;

    M2 = zeros(D, s, D, s) ;
    MN2 = zeros(D, s, D, s) ;


    I= operator_2=kron(kron(kron(u,u),u),u)


    M1[1,:,1,:] = I ; M1[ D,:,D,:] = I
    MN1[1,:,1,:] = I ;    MN1[ D,:,D,:] = I
    M2[1,:,1,:] = I ; M2[ D,:,D,:] = I
    MN2[1,:,1,:] = I ;    MN2[ D,:,D,:] = I

    operator_5=kron(kron(kron(pm,u),pm),pp)
    operator_6=kron(kron(kron(pp,u),u),u)

    operator_7=kron(kron(kron(u,pm),pp),pm)
    operator_8=kron(kron(kron(u,pp),u),u)


    operator_5_dag=kron(kron(kron(pp,u),pp),pm)
    operator_6_dag=kron(kron(kron(pm,u),u),u)

    operator_7_dag=kron(kron(kron(u,pp),pm),pp)
    operator_8_dag=kron(kron(kron(u,pm),u),u)


    # first plaquette
    M1[1,:,2,:] = operator_5
    M1[2,:,D,:] = operator_6

    M1[1,:,3,:] = operator_5_dag
    M1[3,:,D,:] = operator_6_dag

    MN1[2,:,D,:]  = operator_6
    MN1[3,:,D,:] = operator_6_dag

    # second plaquette
    M2[1,:,2,:] = operator_7
    M2[2,:,D,:] = operator_8

    M2[1,:,3,:] = operator_7_dag
    M2[3,:,D,:] = operator_8_dag

    MN2[2,:,D,:] = operator_8
    MN2[3,:,D,:] = operator_8_dag


    mpo= [ M1[1:1,:,:,:] ]

    if sublattice=="A"
        for site = 2:N-1;
            if site % 2 == 0
                push!(mpo, M1[:,:,:,:])
            else
                push!(mpo, M2[:,:,:,:])
            end
        end

        if (N-1)%2==0
            push!(mpo, MN1[:,:,D:D,:])
        else
            push!(mpo, MN2[:,:,D:D,:])
        end
    elseif sublattice=="B"
        for site = 2:N-1;
            if site % 2 == 0
                push!(mpo, M2[:,:,:,:])
            else
                push!(mpo, M1[:,:,:,:])
            end
        end

        if (N-1)%2 == 0
            push!(mpo, MN2[:,:,D:D,:])
        else
            push!(mpo, MN1[:,:,D:D,:])
        end
    else
        print("\nI dont know this wired sublattice called: '", sublattice, "'\n")
        print("Please specify with sublattice =['A' or 'B']")
        return
    end

    # construct MPO


    return mpo
end




#print(operator_Oflip(10))

#print(operator_sublattice_energy(10,"A"))