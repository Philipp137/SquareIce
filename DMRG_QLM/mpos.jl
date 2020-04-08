module mpos
#export plaquette_operators
using LinearAlgebra, TensorOperations, KrylovKit

#################################################
# Paulimatrices and projectors:
sp = [0. 1.; 0. 0.]     # sigma plus
sm = [0. 0.; 1. 0.]     # sigma minus
sz = [1. 0.; 0. -1.]    # sigma z
u = [1. 0.; 0. 1.]      # Identity
pp = [1. 0.; 0. 0.]     # p plus
pm = [0. 0.; 0. 1.]     # p minus
#################################################


function plaquette_operators(Ly::Int=2 , theta=0.)

    if abs(theta)>0
        phi=exp(im*theta)
        phi_star=exp(-1*im*theta)
    else
        phi=phi_star=1
    end

    if Ly == 2
        plaquette_ops=[# first plaquette
                       -1*phi*kron(kron(kron(sm,u),sm),sp),
                       kron(kron(kron(sp,u),u),u),
                       # second plaquette
                       -1*phi*kron(kron(kron(u,sm),sp),sm),
                       kron(kron(kron(u,sp),u),u),
                       ## dagger operators
                       -1*phi_star*kron(kron(kron(sp,u),sp),sm),
                       kron(kron(kron(sm,u),u),u),
                       -1*phi_star*kron(kron(kron(u,sp),sm),sp),
                       kron(kron(kron(u,sm),u),u)]
    elseif Ly == 3
        plaquette_ops=[# first plaquette
                       -1*phi*kron(kron(kron(kron(kron(sm,u),u),sm),sp),u),
                       kron(kron(kron(kron(kron(sp,u),u),u),u),u),
                       # second plaquette
                       -1*phi*kron(kron(kron(kron(kron(u,sm),u),u),sm),sp),
                       kron(kron(kron(kron(kron(u,sp),u),u),u),u),
                       # third plaquette
                       -1*phi*kron(kron(kron(kron(kron(u,u),sm),sp),u),sm),
                       kron(kron(kron(kron(kron(u,u),sp),u),u),u),
                       ## dagger operators
                       -1*phi_star*kron(kron(kron(kron(kron(sp,u),u),sp),sm),u),
                       kron(kron(kron(kron(kron(sm,u),u),u),u),u),
                       -1*phi_star*kron(kron(kron(kron(kron(u,sp),u),u),sp),sm),
                       kron(kron(kron(kron(kron(u,sm),u),u),u),u),
                       -1*phi_star*kron(kron(kron(kron(kron(u,u),sp),sm),u),sp),
                       kron(kron(kron(kron(kron(u,u),sm),u),u),u)]
    else
        error("Ly either 2 or 3")
    end
    return plaquette_ops
end

function interaction_operators(Ly::Int=2, coupling = -1.)

    if Ly == 2
        ops=[#first interaction term
            coupling*kron(kron(kron(pm,u),pm),pp),
            kron(kron(kron(pp,u),u),u),
            # second interaction term
            coupling*kron(kron(kron(u,pm),pp),pm),
            kron(kron(kron(u,pp),u),u),
            ## dagger ops
            coupling*kron(kron(kron(pp,u),pp),pm),
            kron(kron(kron(pm,u),u),u),
            #
            coupling*kron(kron(kron(u,pp),pm),pp),
            kron(kron(kron(u,pm),u),u),
            # projector
            kron(kron(kron(u,u),pp),pm)
            ]
    elseif Ly == 3
        ops=[# first interaction
                coupling*kron(kron(kron(kron(kron(pm,u),u),pm),pp),u),
                kron(kron(kron(kron(kron(pp,u),u),u),u),u),
               # second interaction
                coupling*kron(kron(kron(kron(kron(u,pm),u),u),pm),pp),
                kron(kron(kron(kron(kron(u,pp),u),u),u),u),
                # third interaction
                coupling*kron(kron(kron(kron(kron(u,u),pm),pp),u),pm),
                kron(kron(kron(kron(kron(u,u),pp),u),u),u),
                ## dagger operators
                coupling*kron(kron(kron(kron(kron(pp,u),u),pp),pm),u),
                kron(kron(kron(kron(kron(pm,u),u),u),u),u),
                coupling*kron(kron(kron(kron(kron(u,pp),u),u),pp),pm),
                kron(kron(kron(kron(kron(u,pm),u),u),u),u),
                coupling*kron(kron(kron(kron(kron(u,u),pp),pm),u),pp),
                kron(kron(kron(kron(kron(u,u),pm),u),u),u),
                # projector
                kron(kron(kron(kron(kron(u,u),u),pp),pm),pp)
                ]
    else
        error("Ly either 2 or 3")
    end
    return ops
end


function chemical_potential_operators(Ly::Int=2, mu=-1)

    if Ly == 2
        ops=[
            kron(kron(kron(sz,u),u),u),
            kron(kron(kron(u,sz),u),u),
            kron(kron(kron(u,u),sz),u),
            kron(kron(kron(u,u),u),sz)
            ]
    elseif Ly == 3
        ops=[
            kron(kron(kron(kron(kron(sz,u),u),u),u),u),
            kron(kron(kron(kron(kron(u,sz),u),u),u),u),
            kron(kron(kron(kron(kron(u,u),sz),u),u),u),
            kron(kron(kron(kron(kron(u,u),u),sz),u),u),
            kron(kron(kron(kron(kron(u,u),u),u),sz),u),
            kron(kron(kron(kron(kron(u,u),u),u),u),sz)
            ]
    else
        error("Ly either 2 or 3")
    end
    return mu*ops
end


function mpoqlm(N::Int, Ly::Int=2, coupling=-1. , mu=-1 , theta= 0.)

    Nlinks = Ly * 2
    Nstates = 2^(Nlinks) # number of states for one chain element

    # operator dimension / size
    if Ly == 2
        D = 10
    elseif Ly == 3
        D = 12
    else
            error("Ly either 2 or 3")
    end

    # initialice
    M = zeros(D, Nstates, D, Nstates) ;
    MN = zeros(D, Nstates, D, Nstates) ;

    I= u
    for k = 1 : Nlinks-1
        I= kron(I,u)
    end
    M[1,:,1,:] = I ;    M[ D,:,D,:] = I
    MN[1,:,1,:] = I ;    MN[ D,:,D,:] = I
    M = convert(Array{Complex{Float64},4},M)
    MN = convert(Array{Complex{Float64},4},MN)

    ##########
    # compute operators with respect to local basis
    ##########
    # Plaquette Operators inside magnetic field with angle theta
    plq_ops   = plaquette_operators(Ly , theta)
    # Interaction Terms and projector
    Inter_ops = interaction_operators(Ly, coupling)
    # chemical potential
    mu_ops    = chemical_potential_operators(Ly, mu)

    M[1,:,D,:] = sum(mu_ops[Ly+1:end])

    # Plaquette Terms
    for ip = 1:Nlinks;
        M[1,:,ip+1,:] = plq_ops[2*ip-1]
        M[ip+1,:,D,:] = plq_ops[2*ip]
    end

    # Interaction terms
    for k = 1:Nlinks;
        M[1,:,5+k,:] = Inter_ops[2*k-1]
        M[5+k,:,D,:] = Inter_ops[2*k]
    end

    MN[1,:,D,:] = sum(mu_ops[1:Ly])

    MN[1,:,D,:] = Inter_ops[end] # last element of Inter_ops contains projector

    for ip = 1:Nlinks;
        MN[ip+1,:,D,:] = plq_ops[2*ip]
    end

    # Interaction terms
    for k = 1:Nlinks;
        MN[5+k,:,D,:] = Inter_ops[2*k]
    end

    # Construct Matrix Product Operator
    mpo= [ M[1:1,:,:,:] ]

    for site = 2:N-1;
        push!(mpo, M[:,:,:,:])
    end

    push!(mpo, MN[:,:,D:D,:])

    return mpo
end

end
