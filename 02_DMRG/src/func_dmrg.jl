using LinearAlgebra, TensorOperations, KrylovKit, Printf, DelimitedFiles

randisometry(T, d1, d2) = d1 >= d2 ? Matrix(qr!(randn(T, d1, d2)).Q) : Matrix(lq!(randn(T, d1, d2)).Q)
randisometry(d1, d2) = randisometry(Float64, d1, d2)
randisometry(dims::Dims{2}) = randisometry(dims[1], dims[2])
randisometry(T, dims::Dims{2}) = randisometry(T, dims[1], dims[2])

"""
    randmps(physdims::NTuple{N,Int}, Dmax::Int, [T::Type{<:Number} = Float64])
    randmps(N::Int, d::Int, Dmax::Int, [T::Type{<:Number} = Float64])

Construct a random right canonical MPS for a system with `N`, where site `n` has local Hilbert
space dimension `physdims[n]` (first method) or `d` (second method), and the maximal bond
dimension is `Dmax`. Entries of the MPS tensors will be of type `T`, defaulting to `Float64`.
"""
function randmps(physdims::Dims{N}, Dmax::Int, T::Type{<:Number} = Float64) where {N}
    bonddims = Vector{Int}(undef, N+1)
    bonddims[1] = 1
    bonddims[N+1] = 1
    Nhalf = div(N,2)
    for n = 2:N
        bonddims[n] = min(Dmax, bonddims[n-1]*physdims[n-1])
    end
    for n = N:-1:1
        bonddims[n] = min(bonddims[n], bonddims[n+1]*physdims[n])
    end

    As = Vector{Any}(undef, N)
    for n = 1:N
        d = physdims[n]
        Dl = bonddims[n]
        Dr = bonddims[n+1]
        As[n] = reshape(randisometry(T, Dl, d*Dr), (Dl, d, Dr))
    end
    return As
end
randmps(N::Int, d::Int, Dmax::Int, T = Float64) = randmps(ntuple(n->d, N), Dmax, T)

function applyH1(AC, FL, FR, M)
    @tensor HAC[α,s,β] := FL[α,a,α']*AC[α',s',β']*M[a,s,b,s']*FR[β',b,β]
end

function applyH2(AAC, FL, FR, M1, M2)
    @tensor HAAC[α,s1,s2,β] := FL[α,a,α']*AAC[α',s1',s2',β']*M1[a,s1,b,s1']*M2[b,s2,c,s2']*FR[β',c,β]
end

function svdtrunc(A; truncdim = max(size(A)...), truncerr = 0.)
    F = svd(A)
    d = min(truncdim, count(F.S .>= truncerr))
    return F.U[:,1:d], diagm(0=>F.S[1:d]), F.Vt[1:d, :]
end

function updateleftenv(A, M, FL)
    @tensor FL[α,a,β] := FL[α',a',β']*A[β',s',β]*M[a',s,a,s']*conj(A[α',s,α])
end

function updaterightenv(A, M, FR)
    @tensor FR[α,a,β] := A[α,s',α']*FR[α',a',β']*M[a,s,a',s']*conj(A[β,s,β'])
end


function dmrg1sweep!(A, M, F = nothing; verbose = true, kwargs...)
    N = length(A)

     if F == nothing
         F = Vector{Any}(undef, N+2)
         F[1] = fill!(similar(M[1], (1,1,1)), 1)
         F[N+2] = fill!(similar(M[1], (1,1,1)), 1)
        for k = N:-1:1
            F[k+1] = updaterightenv(A[k], M[k], F[k+2])
        end
    end

    AC = A[1]
    for k = 1:N-1
        Es, ACs, info = eigsolve(x->applyH1(x, F[k], F[k+2], M[k]), AC, 1, :SR; ishermitian = true, kwargs...)
        AC = ACs[1]
        E = Es[1]

        verbose && println("Sweep L2R: site $k -> energy $E")

        AL, C = qr(reshape(AC, size(AC,1)*size(AC,2), :))
        A[k] = reshape(Matrix(AL), size(AC))
        F[k+1] = updateleftenv(A[k], M[k], F[k])

        @tensor AC[-1,-2,-3] := C[-1,1] * A[k+1][1,-2,-3]
    end
    k = N
    Es, ACs, info = eigsolve(x->applyH1(x, F[k], F[k+2], M[k]), AC, 1, :SR; ishermitian = true, kwargs...)
    AC = ACs[1]
    E = Es[1]
    verbose && println("Sweep L2R: site $k -> energy $E")
    for k = N-1:-1:1
        C, AR = lq(reshape(AC, size(AC,1), :))
        # it's actually better to do qr of transpose and transpose back

        A[k+1] = reshape(Matrix(AR), size(AC))
        F[k+2] = updaterightenv(A[k+1], M[k+1], F[k+3])

        @tensor AC[:] := A[k][-1,-2,1] * C[1,-3]
        Es, ACs, info = eigsolve(x->applyH1(x, F[k], F[k+2], M[k]), AC, 1, :SR; ishermitian = true, kwargs...)
        AC = ACs[1]
        E = Es[1]
        verbose && println("Sweep R2L: site $k -> energy $E")
    end
    A[1] = AC
    return E, A, F
end



function dmrg2sweep!(A, M, F = nothing; verbose = true, truncdim = 200, truncerr = 1e-6, kwargs...)
    N = length(A)
    if F == nothing
        F = Vector{Any}(undef, N+2)
        F[1] = fill!(similar(M[1], (1,1,1)), 1)
        F[N+2] = fill!(similar(M[1], (1,1,1)), 1)
        for k = N:-1:1
            F[k+1] = updaterightenv(A[k], M[k], F[k+2])
        end
    end

    AC = A[1]
    for k = 1:N-2
        @tensor AAC[-1,-2,-3,-4] := AC[-1,-2,1]*A[k+1][1,-3,-4]
        Es, AACs, info = eigsolve(x->applyH2(x, F[k], F[k+3], M[k], M[k+1]), AAC, 1, :SR; ishermitian = true, kwargs...)
        AAC = AACs[1]
        E = Es[1]

        verbose && println("Sweep L2R: site $(k:k+1) -> energy $E")

        AL, S, V = svdtrunc(reshape(AAC, size(AAC,1)*size(AAC,2), :); truncdim = truncdim, truncerr = truncerr)
        A[k] = reshape(AL, size(AC, 1), size(AC, 2), :)
        F[k+1] = updateleftenv(A[k], M[k], F[k])

        AC = reshape(S*V, size(S,1), size(A[k+1], 2), size(A[k+1], 3))
    end

    k = N-1
    @tensor AAC[-1,-2,-3,-4] := AC[-1,-2,1]*A[k+1][1,-3,-4]
    Es, AACs, info = eigsolve(x->applyH2(x, F[k], F[k+3], M[k], M[k+1]), AAC, 1, :SR; ishermitian = true, kwargs...)
    AAC = AACs[1]
    E = Es[1]
    verbose && println("Sweep L2R: site $(k:k+1) -> energy $E")

    for k = N-1:-1:2
        U, S, AR = svdtrunc(reshape(AAC, size(AAC,1)*size(AAC,2), :); truncdim = truncdim, truncerr = truncerr)

        A[k+1] = reshape(AR, size(AR, 1), size(AAC, 3), size(AAC, 4))
        F[k+2] = updaterightenv(A[k+1], M[k+1], F[k+3])

        AC = reshape(U*S, size(AAC,1), size(AAC,2), size(S,2))
        @tensor AAC[:] := A[k-1][-1,-2,1] * AC[1,-3,-4]

        Es, AACs, info = eigsolve(x->applyH2(x, F[k-1], F[k+2], M[k-1], M[k]), AAC, 1, :SR; ishermitian = true, kwargs...)
        AAC = AACs[1]
        E = Es[1]
        verbose && println("Sweep R2L: site $(k-1:k) -> energy $E")
    end
    k = 1
    U, S, AR = svdtrunc(reshape(AAC, size(AAC,1)*size(AAC,2), :); truncdim = truncdim, truncerr = truncerr)

    A[k+1] = reshape(AR, size(AR, 1), size(AAC, 3), size(AAC, 4))
    F[k+2] = updaterightenv(A[k+1], M[k+1], F[k+3])

    AC = reshape(U*S, size(AAC,1), size(AAC,2), size(S,2))
    A[1] = AC
    return E, A, F
end



function dmrgconvergence!(A, M , F = nothing ;  verbose = true, kwargs...)
    N = length(A)

    if F == nothing
        F = Vector{Any}(undef, N+2)
        F[1] = fill!(similar(M[1], (1,1,1)), 1)
        F[N+2] = fill!(similar(M[1], (1,1,1)), 1)
        for k = N:-1:1
            F[k+1] = updaterightenv(A[k], M[k], F[k+2])
        end
    end

    max_sweep=1000000
    conv = 1.0e-8
    E = Vector{Float64}(undef, max_sweep)
    counter=2
    E[1], A,  F = dmrg1sweep!( A, M; verbose = false);
#    println("1")
#    println(E[1])
    E[2]=1.

    while  abs(E[counter]-E[counter-1]) > conv
        counter+=1
        E[counter], A, F = dmrg1sweep!(A, M, F; verbose = false);
#        println(counter-1)
#        println(E[counter])
    end
    X=E[counter]
    return X , A, F
end

function measure1siteoperator(A, O)
    N = length(A)
    ρ = ones(eltype(A[1]), 1, 1)
    expval = Vector{Complex{Float64}}(undef, N)
    for k = 1:N
        @tensor v = scalar(ρ[a,b]*A[k][b,s,c]*O[s',s]*conj(A[k][a,s',c]))
        expval[k] = v
        @tensor ρ[a,b] := ρ[a',b']*A[k][b',s,b]*conj(A[k][a',s,a])
    end
    return expval
end

function measure_two_NNN_site_operator(A, O_1, O_2, k_1, k_2)
    N = length(A)
    ρ = ones(eltype(A[1]), 1, 1)

    for k = 1:N
        if k == k_2
            @tensor ρ[a, b] :=
                ρ[a', b'] * A[k][b', s, b] * O_2[s', s] * conj(A[k][a', s', a])
        elseif k == k_1
            @tensor ρ[a, b] :=
                ρ[a', b'] * A[k][b', s, b] * O_1[s', s] * conj(A[k][a', s', a])
        else
            @tensor ρ[a, b] := ρ[a', b'] * A[k][b', s, b] * conj(A[k][a', s, a])
        end
    end
    eta = ones(eltype(A[1]), 1, 1)
    @tensor v = scalar(ρ[a, b] * eta[a, b])

    return v
end


function measure_mpo!(A, M )
    N = length(A)
    F = Vector{Any}(undef, N+2)
    F[1] = fill!(similar(M[1], (1,1,1)), 1)
    F[N+2] = fill!(similar(M[1], (1,1,1)), 1)
    for k = N:-1:1
        F[k+1] = updaterightenv(A[k], M[k], F[k+2])
    end
    for k = 1:N
        F[k+1] = updateleftenv(A[k], M[k], F[k])
    end
    FL=F[N+1]
    FR=F[N+2]
   @tensor E = scalar(  FL[α,a,α']*FR[α',a,α]   )
    return E
    end



function dmrgconvergence_in_D!(s, D, D_max , A, M , F = nothing ;  verbose = true, kwargs...)
    N = length(A)

    A = randmps(N, s, D);

    max_sweep=100
    conv = 1.0e-8
    E = Vector{Float64}(undef, max_sweep)
    counter=2

    E[2], A, F = dmrg1sweep!( A, M; verbose = false);

    B=[]
    G=[]
#    push!(A,B)
#    push!(F,G)

    E[1]=E[2]+1.
    #println("$N  $D    $(E[2])   ")
    while  abs(E[counter]-E[counter-1]) > conv
        D=2*D
        counter+=1
        E[counter], A, F = dmrg2sweep!(A , M ; verbose = false, truncdim = D , truncerr = 1e-10 )
        E[counter], A, F  = dmrgconvergence!(A, M , F  ; verbose = true);
#        push!(A,B)
#        push!(F,G)
    #    println("$N  $D    $(E[counter])   ")
    end
    return E , A, F, counter
#    return E , B, G

end

function measure_correlator!(A, Operator, Lx)
    offset = 3 # offset from the boundaries
    corr =  Vector{Float64}(undef,Lx-2*offset)
    distance =  Vector{Float64}(undef,Lx-2*offset)
    pos_1 = offset
    for r = 1:Lx-2*offset
         pos_2 = pos_1 + r
         distance[r] = r
         corr[r] = measure_two_NNN_site_operator(A, Operator, Operator, pos_1, pos_2)
    end
    return corr, distance
end

function dmrgconvergence_in_D_and_measure_op!(coupling_interaction ,chemical_potential  ,  theta  , s, D, D_max, A, M, F = nothing ;  verbose = false, kwargs...)
    N = length(A)

    A = randmps(N, s, D);

    max_sweep = 100
    conv = 1.0e-8
    E = Vector{Float64}(undef, max_sweep)
    counter = 2
    sp = [0. 1.; 0. 0.]
    sm = [0. 0.; 1. 0.]
    sz = [1. 0.; 0. -1.]
    u = [1. 0.; 0. 1.]
    pp = [1. 0.; 0. 0.]
    pm = [0. 0.; 0. 1.]
    O = kron(kron(kron(sz, sz), u), u)

    winding_number = measure1siteoperator(A, O)
    winding_number = deleteat!(winding_number, N)
    # define operators
    Oflipp_mpo = operator_flipp(N)
    Oflip_mpo = operator_Oflip(N)
    EA_mpo=operator_sublattice_energy(N,"A")
    EB_mpo=operator_sublattice_energy(N,"B")
    MA_mpo = chess_operator_down( N )
    MB_mpo = chess_operator_up( N )
    MB2_mpo = chess_operator_squared( N ) # square of the chess operator <M_B^2> = <M_B^4>
    # measure operators
    EA    = measure_mpo!(A,EA_mpo)
    EB    = measure_mpo!(A,EB_mpo)
    MA    = measure_mpo!(A,MA_mpo)
    MB    = measure_mpo!(A,MB_mpo)
    MB2    = measure_mpo!(A,MB2_mpo) # square of the chess operator <M_B^2> = <M_B^4>
    Oflip = measure_mpo!(A,Oflip_mpo)
    Oflipp = measure_mpo!(A, Oflipp_mpo)
    charge = measure_charge(A)
    entropy = Von_Neumann_entropy(A)

    E[2], A, F  = dmrgconvergence!(A, M, F  ; verbose = true);
    @printf("%4s %6s %8s %5s %5s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n",
            "Lx","lambda","mu_y", "theta", "bondD ","Energy_GS","winding","eA","eB",
            "MA","MB","U4","Oflip","Oflipp", "charge")
    # @printf("%4d %6.2f %8.2f %5.2f %5d %10f %10f %10f %10f %10f %10f %10f %10f \n",(N - 1),
    #         coupling_interaction,chemical_potential,theta,D,(E[counter-1]),(real(sum(winding_number))),
    #         (real(EA)),(real(EB)),(real(MA)),(real(MB)),(real(Oflip)),(real(Oflipp)))
    E[1] = E[2] + 1.
    if verbose
        println("$N  $D    $(E[2])   ")
    end

    #println("Number_Plaquettes\ncoupling\nchemical\ntheta\nBond_dimention       Energy_GS                        winding_number                 flipp ")

    D0 = D
    while  abs(E[counter] - E[counter - 1]) > conv
        D = D + D0
        counter += 1
        E[counter], A, F = dmrg2sweep!(A, M ; verbose = false, truncdim = D , truncerr = 1e-10)
        E[counter], A, F  = dmrgconvergence!(A, M, F  ; verbose = true);

        # Measure
        winding_number = measure1siteoperator(A, O)
        winding_number = deleteat!(winding_number, N)
        EA    = measure_mpo!(A,EA_mpo)
        EB    = measure_mpo!(A,EB_mpo)
        MA    = measure_mpo!(A,MA_mpo)
        MB    = measure_mpo!(A,MB_mpo)
        MB2    = measure_mpo!(A,MB2_mpo)
        U4 = 1 - MB2/(3*MB2^2)
        Oflip = measure_mpo!(A,Oflip_mpo)
        Oflipp = measure_mpo!(A, Oflipp_mpo)
        charge = measure_charge(A)
        entropy = Von_Neumann_entropy(A)
        # print to console
        @printf("%4d %6.2f %8.2f %5.2f %5d %10f %10f %10f %10f %10f %10f %10f %10f %10f %10f\n",(N - 1),
                coupling_interaction,chemical_potential[2],theta,D,(E[counter-1]),
                (real(sum(winding_number))),(real(EA)),(real(EB)),(real(MA)),
                (real(MB)),(real(U4)),(real(Oflip)),(real(Oflipp)),real(charge))
        if verbose
            println("$N  $D    $(E[counter])   ")
        end
    end


    return E, A, F , counter
#    return E , B, G

end



function measure_op!(coupling , mu  ,  theta  , s, D, D_max, A, M, F = nothing ;  verbose = false)
    N = length(A)
    max_sweep = 100
    conv = 1.0e-8
    E = Vector{Float64}(undef, max_sweep)
    counter = 2
    sp = [0. 1.; 0. 0.]
    sm = [0. 0.; 1. 0.]
    sz = [1. 0.; 0. -1.]
    u = [1. 0.; 0. 1.]
    pp = [1. 0.; 0. 0.]
    pm = [0. 0.; 0. 1.]
    O = kron(kron(kron(sz, sz), u), u)
    winding_number = measure1siteoperator(A, O)
    winding_number = deleteat!(winding_number, N)

    # define operators
    Oflipp_mpo = operator_flipp(N)
    Oflip_mpo  = operator_Oflip(N)
    EA_mpo     = operator_sublattice_energy(N,"A")
    EB_mpo     = operator_sublattice_energy(N,"B")
    MA_mpo     = chess_operator_down( N )
    MB_mpo     = chess_operator_up( N )
    MB2_mpo    = chess_operator_squared( N ) # square of the chess operator <M_B^2> = <M_B^4>
    # measure operators
    EA         = measure_mpo!(A,EA_mpo)
    EB         = measure_mpo!(A,EB_mpo)
    MA         = measure_mpo!(A,MA_mpo)
    MB         = measure_mpo!(A,MB_mpo)
    MB2        = measure_mpo!(A,MB2_mpo) # square of the chess operator <M_B^2> = <M_B^4>
    Oflip      = measure_mpo!(A,Oflip_mpo)
    Oflipp     = measure_mpo!(A, Oflipp_mpo)
    charge     = measure_charge(A)
    entropy    = Von_Neumann_entropy(A)
    U4 = 1 - MB2/(3*MB2^2)
    E[2], A, F = dmrg1sweep!( A, M; verbose = false);
    @printf("%4s %6s %8s %5s %5s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n",
            "Lx","lambda","mu_y", "theta", "bondD ","Energy_GS","winding","eA","eB",
            "MA","MB","U4","Oflip","Oflipp", "charge")
    @printf("%4d %6.2f %8.2f %5.2f %5d %10f %10f %10f %10f %10f %10f %10f %10f %10f %10f\n",(N - 1),
                    coupling, mu[2],theta,D,(E[2]),
                    (real(sum(winding_number))),(real(EA)),(real(EB)),(real(MA)),
                    (real(MB)),(real(U4)),(real(Oflip)),(real(Oflipp)),real(charge))

    correlator ,distance = measure_correlator!(A, O, N-1)
    fname = string("winding_correlator.txt" )
    @printf("Saving Corr-fun to: %s\n", fname)
    open(fname, "w") do io
        writedlm(io, [distance correlator])
    end
end
