
using LinearAlgebra, TensorOperations, KrylovKit

function exp_value(A, O)
    @tensor v = scalar(A[a]*O[a,b]*conj(A[b] ))
    return v
end

function subspace_projector( )

    sz = [1. 0.; 0. -1.]
    u = [1. 0.; 0. 1.]
    I= operator_2=kron(kron(kron(kron(kron(u,u),u),u),u),u)
    Sz_1= operator_2=kron(kron(kron(kron(kron(sz,u),u),u),u),u)
    Sz_2= operator_2=kron(kron(kron(kron(kron(u,sz),u),u),u),u)
    Sz_3= operator_2=kron(kron(kron(kron(kron(u,u),sz),u),u),u)
    Sz_4= operator_2=kron(kron(kron(kron(kron(u,u),u),sz),u),u)
    Sz_5= operator_2=kron(kron(kron(kron(kron(u,u),u),u),sz),u)
    Sz_6= operator_2=kron(kron(kron(kron(kron(u,u),u),u),u),sz)
    
    GL=Sz_1 + - Sz_4 - Sz_3 + Sz_5

    GR= Sz_2 + - Sz_4 - Sz_3 + Sz_6

    
    GTOT= GL*GL + GR*GR
    F =    eigen(GTOT) 

    M3_test = zeros(64,64)

    i=1
    j=0
    while i < length(F.values) +1
        if abs( F.values[i] ) < 0.0000001  &&  abs(exp_value(F.vectors[:, i],GR) ) < 0.0000001 &&  abs(exp_value(F.vectors[:, i],GL) ) < 0.0000001 &&  abs(exp_value(F.vectors[:, i],I) ) > 0.0000001
            j = j+1
            M3_test[:,j]=F.vectors[:, i]
#            #println(M3_test[:,j])
        end

        i = i+1
    end

    M3=zeros(64,j)
    y=1
    while y < j+1
        M3[:,y]=M3_test[:,y]
#        #println(M3[:,y])
        y=y+1
    end

    return  M3
end


function test_subspace_projected_operators( )

    projector=subspace_projector()
    y=1
    while y < 11
        #println(projector[:,y])
        y=y+1
    end

    s=10
    D = 14
    sp = [0.  1. ;  0.  0.]
    sx = [0.  1. ;  0.  0.]
    sm = [0.  0. ;  1.  0.]
    sz = [1.  0. ;  0. -1.]
    u =  [1.  0. ;  0.  1.]

    I= kron(kron(kron(kron(kron(u,u),u),u),u),u)
    @tensor I[α,β] := projector[α',α]*I[α',β']*projector[β',β]

    #println(I)

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
    
    
    Sx_6=  kron(kron(kron(kron(kron(u,u),u),u),u),sx)

    projector_1=kron(kron(kron(kron(kron((u+sz)*0.5,(u-sz)*0.5),u),u),u),u)
    projector_N=kron(kron(kron(kron(kron(u,u),u),u),(u+sz)*0.5),(u-sz)*0.5)
    GL=Sz_1 + - Sz_4 - Sz_3 + Sz_5

    GR= Sz_2 + - Sz_4 - Sz_3 + Sz_6

    @tensor Sz_1_Pjrojected[α,β] := projector[α',α]*Sz_1[α',β']*projector[β',β]
    
    #println(Sz_1_Pjrojected)
    GTOT= GL*GL + GR*GR

    operator_1=Sp_5*Sm_6*Sm_3
    @tensor operator_1_Pjrojected[α,β] := projector[α',α]*operator_1[α',β']*projector[β',β]
    println(operator_1_Pjrojected)
    
    GTOT= GL*GL + GR*GR


end


test_subspace_projected_operators()



#println(subspace_projector())