using LinearAlgebra, TensorOperations, KrylovKit

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


function updateleftenv(A, M, FL)
    @tensor FL[α,a,β] := FL[α',a',β']*A[β',s',β]*M[a',s,a,s']*conj(A[α',s,α])
end

function updaterightenv(A, M, FR)
    @tensor FR[α,a,β] := A[α,s',α']*FR[α',a',β']*M[a,s,a',s']*conj(A[β,s,β'])
end


function dmrg1sweep!(A, H, F = nothing; verbose = true, kwargs...)
    N = length(A)
    
    if F == nothing
        F = Vector{Any}(undef, N+2)
        F[1] = fill!(similar(H[1], (1,1,1)), 1)
        F[N+2] = fill!(similar(H[1], (1,1,1)), 1)
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


function dmrgconvergence!(A, M ;  verbose = true, kwargs...)
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

