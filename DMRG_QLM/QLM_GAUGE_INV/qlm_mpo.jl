
using LinearAlgebra, TensorOperations, KrylovKit

include("diagonalize_gauss_operator.jl")

function mpoqlm_6_link_base(N::Int;   coupling=1. )
    sp = [0. 1.; 0. 0.]
    sm = [0. 0.; 1. 0.]
    sz = [1. 0.; 0. -1.]
    u = [1. 0.; 0. 1.]


    s=64
    D = 14
    M = zeros(D, s, D, s) ;

    I= operator_2=kron(kron(kron(kron(kron(u,u),u),u),u),u)


    M[1,:,1,:] = I ;    M[ D,:,D,:] = I

    Sz_1= kron(kron(kron(kron(kron(sz,u),u),u),u),u)
    Sz_2= kron(kron(kron(kron(kron(u,sz),u),u),u),u)
    Sz_3= kron(kron(kron(kron(kron(u,u),sz),u),u),u)
    Sz_4= kron(kron(kron(kron(kron(u,u),u),sz),u),u)
    Sz_5= kron(kron(kron(kron(kron(u,u),u),u),sz),u)
    Sz_6= kron(kron(kron(kron(kron(u,u),u),u),u),sz)

    Sp_1= operator_2=kron(kron(kron(kron(kron(sp,u),u),u),u),u)
    Sp_2= operator_2=kron(kron(kron(kron(kron(u,sp),u),u),u),u)
    Sp_3= operator_2=kron(kron(kron(kron(kron(u,u),sp),u),u),u)
    Sp_4= operator_2=kron(kron(kron(kron(kron(u,u),u),sp),u),u)
    Sp_5= operator_2=kron(kron(kron(kron(kron(u,u),u),u),sp),u)
    Sp_6= operator_2=kron(kron(kron(kron(kron(u,u),u),u),u),sp)

    Sm_1= operator_2=kron(kron(kron(kron(kron(sm,u),u),u),u),u)
    Sm_2= operator_2=kron(kron(kron(kron(kron(u,sm),u),u),u),u)
    Sm_3= operator_2=kron(kron(kron(kron(kron(u,u),sm),u),u),u)
    Sm_4= operator_2=kron(kron(kron(kron(kron(u,u),u),sm),u),u)
    Sm_5= operator_2=kron(kron(kron(kron(kron(u,u),u),u),sm),u)
    Sm_6= operator_2=kron(kron(kron(kron(kron(u,u),u),u),u),sm)

    operator_1=Sp_5*Sm_6*Sm_3
    operator_2=Sp_3*Sm_2*Sp_1

    operator_3=Sm_5*Sp_6*Sm_4
    operator_4=Sm_1*Sp_2*Sp_4

    operator_5=Sp_5*Sm_6*Sm_3*Sz_5
    operator_6=Sp_3*Sm_2*Sp_1*Sz_1

    operator_7=Sm_5*Sp_6*Sm_4*Sz_5
    operator_8=Sm_1*Sp_2*Sp_4*Sz_1

    operator_9=Sp_5*Sm_6*Sm_3*Sz_6
    operator_10=Sp_3*Sm_2*Sp_1*Sz_2

    operator_11=Sm_5*Sp_6*Sm_4*Sz_6
    operator_12=Sm_1*Sp_2*Sp_4*Sz_2

    operator_13=Sp_5*Sm_6*Sm_3*Sz_6*Sz_5
    operator_14=Sp_3*Sm_2*Sp_1*Sz_2*Sz_1

    operator_15=Sm_5*Sp_6*Sm_4*Sz_6*Sz_5
    operator_16=Sm_1*Sp_2*Sp_4*Sz_2*Sz_1

    operator_1_dag=Sm_5*Sp_6*Sp_3
    operator_2_dag=Sm_3*Sp_2*Sm_1

    operator_3_dag=Sm_5*Sp_6*Sm_4
    operator_4_dag=Sm_1*Sp_2*Sp_4

    operator_5_dag=Sm_5*Sp_6*Sp_3*Sz_5
    operator_6_dag=Sm_3*Sp_2*Sm_1*Sz_1

    operator_7_dag=Sm_5*Sp_6*Sm_4*Sz_5
    operator_8_dag=Sm_1*Sp_2*Sp_4*Sz_1

    operator_9_dag=Sm_5*Sp_6*Sp_3*Sz_6
    operator_10_dag=Sm_3*Sp_2*Sm_1*Sz_2

    operator_11_dag=Sm_5*Sp_6*Sm_4*Sz_6
    operator_12_dag=Sm_1*Sp_2*Sp_4*Sz_2

    operator_13_dag=Sm_5*Sp_6*Sp_3*Sz_6*Sz_5
    operator_14_dag=Sm_3*Sp_2*Sm_1*Sz_2*Sz_1

    operator_15_dag=Sm_5*Sp_6*Sm_4*Sz_6*Sz_5
    operator_16_dag=Sm_1*Sp_2*Sp_4*Sz_2*Sz_1

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

    M[1,:,10,:] = operator_9
    M[10,:,D,:] = operator_10

    M[1,:,11,:] = operator_11
    M[11,:,D,:] = operator_12

    M[1,:,12,:] = operator_9_dag
    M[12,:,D,:] = operator_10_dag

    M[1,:,13,:] = operator_11_dag
    M[13,:,D,:] = operator_12_dag

    mpo= [ M[1:1,:,:,:] ]

    for site = 2:N-1;
        push!(mpo, M[:,:,:,:])
    end

    push!(mpo, M[:,:,D:D,:])

    return mpo
end


function mpoqlm_gauge_inv(N::Int;   coupling=1. )
    sp = [0.  1. ;  0.  0.]
    sm = [0.  0. ;  1.  0.]
    sz = [1.  0. ;  0. -1.]
    u =  [1.  0. ;  0.  1.]

    projector=subspace_projector()
    s=10
    D = 14
#    D = 6
    M = zeros(D, s, D, s) ;

    I= kron(kron(kron(kron(kron(u,u),u),u),u),u)
    @tensor I[α,β] := projector[α',α]*I[α',β']*projector[β',β]



    M[1,:,1,:] = I ;    M[ D,:,D,:] = I

    Sz_1= kron(kron(kron(kron(kron(sz,u),u),u),u),u)
    Sz_2= kron(kron(kron(kron(kron(u,sz),u),u),u),u)
    Sz_3= kron(kron(kron(kron(kron(u,u),sz),u),u),u)
    Sz_4= kron(kron(kron(kron(kron(u,u),u),sz),u),u)
    Sz_5= kron(kron(kron(kron(kron(u,u),u),u),sz),u)
    Sz_6= kron(kron(kron(kron(kron(u,u),u),u),u),sz)

    Sp_1=  kron(kron(kron(kron(kron(sp,u),u),u),u),u)
    Sp_2=  kron(kron(kron(kron(kron(u,sp),u),u),u),u)
    Sp_3=  kron(kron(kron(kron(kron(u,u),sp),u),u),u)
    Sp_4=  kron(kron(kron(kron(kron(u,u),u),sp),u),u)
    Sp_5=  kron(kron(kron(kron(kron(u,u),u),u),sp),u)
    Sp_6=  kron(kron(kron(kron(kron(u,u),u),u),u),sp)

    Sm_1=  kron(kron(kron(kron(kron(sm,u),u),u),u),u)
    Sm_2=  kron(kron(kron(kron(kron(u,sm),u),u),u),u)
    Sm_3=  kron(kron(kron(kron(kron(u,u),sm),u),u),u)
    Sm_4=  kron(kron(kron(kron(kron(u,u),u),sm),u),u)
    Sm_5=  kron(kron(kron(kron(kron(u,u),u),u),sm),u)
    Sm_6=  kron(kron(kron(kron(kron(u,u),u),u),u),sm)

    projector_1=kron(kron(kron(kron(kron((u+sz)*0.5,(u-sz)*0.5),u),u),u),u)
    projector_N=kron(kron(kron(kron(kron(u,u),u),u),(u+sz)*0.5),(u-sz)*0.5)

    @tensor Sz_1[α,β] := projector[α',α]*Sz_1[α',β']*projector[β',β]
    @tensor Sz_2[α,β] := projector[α',α]*Sz_2[α',β']*projector[β',β]
    @tensor Sz_3[α,β] := projector[α',α]*Sz_3[α',β']*projector[β',β]
    @tensor Sz_4[α,β] := projector[α',α]*Sz_4[α',β']*projector[β',β]
    @tensor Sz_5[α,β] := projector[α',α]*Sz_5[α',β']*projector[β',β]
    @tensor Sz_6[α,β] := projector[α',α]*Sz_6[α',β']*projector[β',β]

    @tensor Sp_1[α,β] := projector[α',α]*Sp_1[α',β']*projector[β',β]
    @tensor Sp_2[α,β] := projector[α',α]*Sp_2[α',β']*projector[β',β]
    @tensor Sp_3[α,β] := projector[α',α]*Sp_3[α',β']*projector[β',β]
    @tensor Sp_4[α,β] := projector[α',α]*Sp_4[α',β']*projector[β',β]
    @tensor Sp_5[α,β] := projector[α',α]*Sp_5[α',β']*projector[β',β]
    @tensor Sp_6[α,β] := projector[α',α]*Sp_6[α',β']*projector[β',β]
    
    @tensor Sm_1[α,β] := projector[α',α]*Sm_1[α',β']*projector[β',β]
    @tensor Sm_2[α,β] := projector[α',α]*Sm_2[α',β']*projector[β',β]
    @tensor Sm_3[α,β] := projector[α',α]*Sm_3[α',β']*projector[β',β]
    @tensor Sm_4[α,β] := projector[α',α]*Sm_4[α',β']*projector[β',β]
    @tensor Sm_5[α,β] := projector[α',α]*Sm_5[α',β']*projector[β',β]
    @tensor Sm_6[α,β] := projector[α',α]*Sm_6[α',β']*projector[β',β]

    @tensor projector_1[α,β] := projector[α',α]*projector_1[α',β']*projector[β',β]
    @tensor projector_N[α,β] := projector[α',α]*projector_N[α',β']*projector[β',β]

    operator_1=Sp_5*Sm_6*Sm_3
    operator_2=Sp_3*Sm_2*Sp_1

    operator_3=Sm_5*Sp_6*Sm_4
    operator_4=Sm_1*Sp_2*Sp_4

    operator_5=Sp_5*Sm_6*Sm_3*Sz_5
    operator_6=Sp_3*Sm_2*Sp_1*Sz_1

    operator_7=Sm_5*Sp_6*Sm_4*Sz_5
    operator_8=Sm_1*Sp_2*Sp_4*Sz_1

    operator_9=Sp_5*Sm_6*Sm_3*Sz_6
    operator_10=Sp_3*Sm_2*Sp_1*Sz_2

    operator_11=Sm_5*Sp_6*Sm_4*Sz_6
    operator_12=Sm_1*Sp_2*Sp_4*Sz_2

    operator_13=Sp_5*Sm_6*Sm_3*Sz_6*Sz_5
    operator_14=Sp_3*Sm_2*Sp_1*Sz_2*Sz_1

    operator_15=Sm_5*Sp_6*Sm_4*Sz_6*Sz_5
    operator_16=Sm_1*Sp_2*Sp_4*Sz_2*Sz_1

    operator_1_dag=Sm_5*Sp_6*Sp_3
    operator_2_dag=Sm_3*Sp_2*Sm_1

    operator_3_dag=Sp_5*Sm_6*Sp_4
    operator_4_dag=Sp_1*Sm_2*Sm_4

    operator_5_dag=Sm_5*Sp_6*Sp_3*Sz_5
    operator_6_dag=Sm_3*Sp_2*Sm_1*Sz_1

    operator_7_dag=Sm_5*Sp_6*Sm_4*Sz_5
    operator_8_dag=Sp_1*Sm_2*Sm_4*Sz_1

    operator_9_dag=Sm_5*Sp_6*Sp_3*Sz_6
    operator_10_dag=Sm_3*Sp_2*Sm_1*Sz_2

    operator_11_dag=Sm_5*Sp_6*Sm_4*Sz_6
    operator_12_dag=Sp_1*Sm_2*Sm_4*Sz_2

    operator_13_dag=Sm_5*Sp_6*Sp_3*Sz_6*Sz_5
    operator_14_dag=Sm_3*Sp_2*Sm_1*Sz_2*Sz_1

    operator_15_dag=Sm_5*Sp_6*Sm_4*Sz_6*Sz_5
    operator_16_dag=Sp_1*Sm_2*Sm_4*Sz_2*Sz_1

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

    M[1,:,10,:] = operator_9
    M[10,:,D,:] = operator_10

    M[1,:,11,:] = operator_11
    M[11,:,D,:] = operator_12

    M[1,:,12,:] = operator_9_dag
    M[12,:,D,:] = operator_10_dag

    M[1,:,13,:] = operator_11_dag
    M[13,:,D,:] = operator_12_dag


    M1=M
    MN=M
    for i = 2:D-1;
        M1[1,:,i,:]=M1[1,:,i,:]*projector_1
        M1[i,:,D,:]=M1[1,:,i,:]*projector_1
    end

    for i = 2:D-1;
        MN[1,:,i,:]=MN[1,:,i,:]*projector_N
        MN[i,:,D,:]=MN[1,:,i,:]*projector_N
    end

    mpo= [ M1[1:1,:,:,:] ]

    for site = 2:N-1;
        push!(mpo, M[:,:,:,:])
    end

    push!(mpo, MN[:,:,D:D,:])

    return mpo
end



function mpoqlm_gauge_inv_2(N::Int;   coupling=1. )
    sp = [0.  1. ;  0.  0.]
    sm = [0.  0. ;  1.  0.]
    sz = [1.  0. ;  0. -1.]
    u =  [1.  0. ;  0.  1.]

    projector=subspace_projector()
    print(size(projector))
    s=10
    D = 14
#    D = 6
    M = zeros(D, s, D, s) ;

    I= kron(kron(kron(kron(kron(u,u),u),u),u),u)
    @tensor I[α,β] := projector[α',α]*I[α',β']*projector[β',β]




    Sz_1= kron(kron(kron(kron(kron(sz,u),u),u),u),u)
    Sz_2= kron(kron(kron(kron(kron(u,sz),u),u),u),u)
    Sz_3= kron(kron(kron(kron(kron(u,u),sz),u),u),u)
    Sz_4= kron(kron(kron(kron(kron(u,u),u),sz),u),u)
    Sz_5= kron(kron(kron(kron(kron(u,u),u),u),sz),u)
    Sz_6= kron(kron(kron(kron(kron(u,u),u),u),u),sz)

    Sp_1=  kron(kron(kron(kron(kron(sp,u),u),u),u),u)
    Sp_2=  kron(kron(kron(kron(kron(u,sp),u),u),u),u)
    Sp_3=  kron(kron(kron(kron(kron(u,u),sp),u),u),u)
    Sp_4=  kron(kron(kron(kron(kron(u,u),u),sp),u),u)
    Sp_5=  kron(kron(kron(kron(kron(u,u),u),u),sp),u)
    Sp_6=  kron(kron(kron(kron(kron(u,u),u),u),u),sp)

    Sm_1=  kron(kron(kron(kron(kron(sm,u),u),u),u),u)
    Sm_2=  kron(kron(kron(kron(kron(u,sm),u),u),u),u)
    Sm_3=  kron(kron(kron(kron(kron(u,u),sm),u),u),u)
    Sm_4=  kron(kron(kron(kron(kron(u,u),u),sm),u),u)
    Sm_5=  kron(kron(kron(kron(kron(u,u),u),u),sm),u)
    Sm_6=  kron(kron(kron(kron(kron(u,u),u),u),u),sm)

  

    operator_1=Sp_5*Sm_6*Sm_3
    @tensor operator_1[α,β] := projector[α',α]*operator_1[α',β']*projector[β',β]
    operator_2=Sp_3*Sm_2*Sp_1
    @tensor operator_2[α,β] := projector[α',α]*operator_2[α',β']*projector[β',β]

    operator_3=Sm_5*Sp_6*Sm_4
    @tensor operator_3[α,β] := projector[α',α]*operator_3[α',β']*projector[β',β]
    operator_4=Sm_1*Sp_2*Sp_4
    @tensor operator_4[α,β] := projector[α',α]*operator_4[α',β']*projector[β',β]

    operator_5=Sp_5*Sm_6*Sm_3*Sz_5
    @tensor operator_5[α,β] := projector[α',α]*operator_5[α',β']*projector[β',β]
    operator_6=Sp_3*Sm_2*Sp_1*Sz_1
    @tensor operator_6[α,β] := projector[α',α]*operator_6[α',β']*projector[β',β]

    operator_7=Sm_5*Sp_6*Sm_4*Sz_5
    @tensor operator_7[α,β] := projector[α',α]*operator_7[α',β']*projector[β',β]
    operator_8=Sm_1*Sp_2*Sp_4*Sz_1
    @tensor operator_8[α,β] := projector[α',α]*operator_8[α',β']*projector[β',β]

    operator_9=Sp_5*Sm_6*Sm_3*Sz_6
    @tensor operator_9[α,β] := projector[α',α]*operator_9[α',β']*projector[β',β]
    operator_10=Sp_3*Sm_2*Sp_1*Sz_2
    @tensor operator_10[α,β] := projector[α',α]*operator_10[α',β']*projector[β',β]

    operator_11=Sm_5*Sp_6*Sm_4*Sz_6
    @tensor operator_11[α,β] := projector[α',α]*operator_11[α',β']*projector[β',β]
    operator_12=Sm_1*Sp_2*Sp_4*Sz_2
    @tensor operator_12[α,β] := projector[α',α]*operator_12[α',β']*projector[β',β]

    operator_13=Sp_5*Sm_6*Sm_3*Sz_6*Sz_5
    @tensor operator_13[α,β] := projector[α',α]*operator_13[α',β']*projector[β',β]
    operator_14=Sp_3*Sm_2*Sp_1*Sz_2*Sz_1
    @tensor operator_14[α,β] := projector[α',α]*operator_14[α',β']*projector[β',β]

    operator_15=Sm_5*Sp_6*Sm_4*Sz_6*Sz_5
    @tensor operator_15[α,β] := projector[α',α]*operator_15[α',β']*projector[β',β]
    operator_16=Sm_1*Sp_2*Sp_4*Sz_2*Sz_1
    @tensor operator_16[α,β] := projector[α',α]*operator_16[α',β']*projector[β',β]

    operator_1_dag=Sm_5*Sp_6*Sp_3
    @tensor operator_1_dag[α,β] := projector[α',α]*operator_1_dag[α',β']*projector[β',β]

    operator_2_dag=Sm_3*Sp_2*Sm_1
    @tensor operator_2_dag[α,β] := projector[α',α]*operator_2_dag[α',β']*projector[β',β]

    operator_3_dag=Sp_5*Sm_6*Sp_4
    @tensor operator_3_dag[α,β] := projector[α',α]*operator_3_dag[α',β']*projector[β',β]
    operator_4_dag=Sp_1*Sm_2*Sm_4
    @tensor operator_4_dag[α,β] := projector[α',α]*operator_4_dag[α',β']*projector[β',β]

    operator_5_dag=Sm_5*Sp_6*Sp_3*Sz_5
    @tensor operator_5_dag[α,β] := projector[α',α]*operator_5_dag[α',β']*projector[β',β]
    operator_6_dag=Sm_3*Sp_2*Sm_1*Sz_1
    @tensor operator_6_dag[α,β] := projector[α',α]*operator_6_dag[α',β']*projector[β',β]

    operator_7_dag=Sm_5*Sp_6*Sm_4*Sz_5
    @tensor operator_7_dag[α,β] := projector[α',α]*operator_7_dag[α',β']*projector[β',β]
    operator_8_dag=Sp_1*Sm_2*Sm_4*Sz_1
    @tensor operator_8_dag[α,β] := projector[α',α]*operator_8_dag[α',β']*projector[β',β]

    operator_9_dag=Sm_5*Sp_6*Sp_3*Sz_6
    @tensor operator_9_dag[α,β] := projector[α',α]*operator_9_dag[α',β']*projector[β',β]
    operator_10_dag=Sm_3*Sp_2*Sm_1*Sz_2
    @tensor operator_10_dag[α,β] := projector[α',α]*operator_10_dag[α',β']*projector[β',β]

    operator_11_dag=Sm_5*Sp_6*Sm_4*Sz_6
    @tensor operator_11_dag[α,β] := projector[α',α]*operator_11_dag[α',β']*projector[β',β]
    operator_12_dag=Sp_1*Sm_2*Sm_4*Sz_2
    @tensor operator_12_dag[α,β] := projector[α',α]*operator_12_dag[α',β']*projector[β',β]

    operator_13_dag=Sm_5*Sp_6*Sp_3*Sz_6*Sz_5
    @tensor operator_13_dag[α,β] := projector[α',α]*operator_13_dag[α',β']*projector[β',β]
    operator_14_dag=Sm_3*Sp_2*Sm_1*Sz_2*Sz_1
    @tensor operator_14_dag[α,β] := projector[α',α]*operator_14_dag[α',β']*projector[β',β]

    operator_15_dag=Sm_5*Sp_6*Sm_4*Sz_6*Sz_5
    @tensor operator_15_dag[α,β] := projector[α',α]*operator_15_dag[α',β']*projector[β',β]
    operator_16_dag=Sp_1*Sm_2*Sm_4*Sz_2*Sz_1
    @tensor operator_16_dag[α,β] := projector[α',α]*operator_16_dag[α',β']*projector[β',β]
    
    projector_1=kron(kron(kron(kron(kron((u+sz)*0.5,(u-sz)*0.5),u),u),u),u)
    projector_N=kron(kron(kron(kron(kron(u,u),u),u),(u+sz)*0.5),(u-sz)*0.5)
    @tensor projector_1[α,β] := projector[α',α]*projector_1[α',β']*projector[β',β]
    @tensor projector_N[α,β] := projector[α',α]*projector_N[α',β']*projector[β',β]
    
    M[1,:,1,:] = I ;    M[ D,:,D,:] = I

    M[1,:,2,:] = operator_1
    println("operator_1")
    println(operator_1)

    M[2,:,D,:] = operator_2

    println("operator_2")
    println(operator_2)

    M[1,:,3,:] = operator_3
    println("operator_3")
    println(operator_3)
    M[3,:,D,:] = operator_4
    println("operator_4")
    println(operator_4)

    M[1,:,4,:] = operator_1_dag
    println("operator_1_dag")
    println(operator_1_dag)
    M[4,:,D,:] = operator_2_dag
    println("operator_2_dag")
    println(operator_2_dag)

    M[1,:,5,:] = operator_3_dag
    println("operator_3_dag")
    println(operator_3_dag)
    M[5,:,D,:] = operator_4_dag
    println("operator_4_dag")
    println(operator_4_dag)

    M[1,:,6,:] = operator_5
    M[6,:,D,:] = operator_6

    M[1,:,7,:] = operator_7
    M[7,:,D,:] = operator_8

    M[1,:,8,:] = operator_5_dag
    M[8,:,D,:] = operator_6_dag

    M[1,:,9,:] = operator_7_dag
    M[9,:,D,:] = operator_8_dag

    M[1,:,10,:] = operator_9
    M[10,:,D,:] = operator_10

    M[1,:,11,:] = operator_11
    M[11,:,D,:] = operator_12

    M[1,:,12,:] = operator_9_dag
    M[12,:,D,:] = operator_10_dag

    M[1,:,13,:] = operator_11_dag
    M[13,:,D,:] = operator_12_dag


    M1=M
    MN=M
    for i = 2:D-1;
        M1[1,:,i,:]=M1[1,:,i,:]*projector_1
        M1[i,:,D,:]=M1[1,:,i,:]*projector_1
    end

    for i = 2:D-1;
        MN[1,:,i,:]=MN[1,:,i,:]*projector_N
        MN[i,:,D,:]=MN[1,:,i,:]*projector_N
    end

    mpo= [ M1[1:1,:,:,:] ]

    
    for site = 2:N-1;
        push!(mpo, M[:,:,:,:])
    end

    push!(mpo, MN[:,:,D:D,:])

    return mpo
end

