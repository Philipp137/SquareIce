using ITensors
using DelimitedFiles

function ITensors.op!(Op::ITensor,
                      ::OpName"Pp",
                      ::SiteType"S=1/2",
                      s::Index)
  Op[s'=>1,s=>1] = 1.

end

function ITensors.op!(Op::ITensor,
                      ::OpName"Pm",
                      ::SiteType"S=1/2",
                      s::Index)
  Op[s'=>2,s=>2] = 1.

end

mutable struct DemoObserver <: AbstractObserver
   energy_tol::Float64
   last_energy::Float64

   DemoObserver(energy_tol=0.0) = new(energy_tol,1000.0)
end

function ITensors.checkdone!(o::DemoObserver;kwargs...)
  sw = kwargs[:sweep]
  energy = kwargs[:energy]
  if abs(energy-o.last_energy)/abs(energy) < o.energy_tol
    println("Stopping DMRG after sweep $sw")
    return true
  end
  # Otherwise, update last_energy and keep going
  o.last_energy = energy
  return false
end


function Hamitonian_qlm(sites, N_x,lambda,J,mu)
    N = N_x*4+2

    ampo = OpSum()

    for j=1:4:N-4
        ampo += lambda,"S+",j+2,"S+",j+4,"S-",j,"S-",j+3
        ampo += lambda,"S+",j+3,"S+",j+5,"S-",j+1,"S-",j+2
        ampo += lambda,"S-",j+2,"S-",j+4,"S+",j,"S+",j+3
        ampo += lambda,"S-",j+3,"S-",j+5,"S+",j+1,"S+",j+2
    end

    for j=1:4:N-4
        ampo += J,"Pm",j+2,"Pm",j+4,"Pp",j,"Pp",j+3
        ampo += J,"Pp",j+2,"Pp",j+4,"Pm",j,"Pm",j+3
        ampo += J,"Pp",j+3,"Pp",j+5,"Pm",j+1,"Pm",j+2
        ampo += J,"Pm",j+3,"Pm",j+5,"Pp",j+1,"Pp",j+2
    end

    for j=1:4:N-1
        ampo += -2*mu,"Sz",j
        ampo += -2*mu,"Sz",j+1
    end

    H = MPO(ampo,sites);
    return H
end;

function Von_N_ent(psi,b)
    orthogonalize!(psi, b)
    U,S,V = svd(psi[b], (linkind(psi, b-1), siteind(psi,b)))
    SvN = 0.0
    for n=1:dim(S, 1)
        p = S[n,n]^2
        SvN -= p * log(p)
    end
    return SvN
end;


function flippability_bot(psi,sites)
    N_x=Int((size(psi)[1]-2)/4)
    flipp_bot=zeros(N_x)
    for j=1:1:N_x
        ampo = OpSum()
        ampo += "Pm",j+2,"Pm",j+4,"Pp",j,"Pp",j+3
        ampo += "Pp",j+2,"Pp",j+4,"Pm",j,"Pm",j+3
        flipp = MPO(ampo,sites);
        flipp_bot[j]=inner(psi' ,flipp,psi)
    end
    return flipp_bot
end

function flippability_up(psi,sites)
    N_x=Int((size(psi)[1]-2)/4)
    flipp_up=zeros(N_x)
    for j=1:1:N_x
        ampo = OpSum()
        ampo += "Pp",j+3,"Pp",j+5,"Pm",j+1,"Pm",j+2
        ampo += "Pm",j+3,"Pm",j+5,"Pp",j+1,"Pp",j+2
        flipp = MPO(ampo,sites);
        flipp_up[j]=inner(psi' ,flipp,psi)
    end
    return flipp_up
end


function dostuff(N_x ,mu; J=-1.0,    lambda=-1. ,storage="./data/")

    N = N_x*4+2
    sites = siteinds("S=1/2",N);
    H=Hamitonian_qlm(sites, N_x,lambda,J,mu)
    etol = 1E-8
    obs = DemoObserver(etol)
    sweeps = Sweeps(200) # number of sweeps is 5
    maxdim!(sweeps,10,40,100,100,200,400,1000) # gradually increase states kept
    cutoff!(sweeps,1E-12) # desired truncation error
    psi0 = randomMPS(sites;linkdims = 20)
    energy, psi = dmrg(H,psi0,sweeps; observer=obs,outputlevel=0);

    magz = expect(psi,"Sz");
    mag=zeros(N_x+1)
    for j=1:Int(N_x+1)
        mag[j] += magz[Int(4*j)-3]+magz[Int(4*j)-2]
    end



    flipp_up=flippability_up(psi,sites)
    flipp_bot=flippability_bot(psi,sites)

    push!(flipp_up,0.);
    push!(flipp_bot,0.);
    E_bot=zeros(size(mag))
    E_up=zeros(size(mag))


    for j=1:Int(N_x)
        E_up[j] += magz[Int(4*j)-1]
    end;

    for j=1:Int(N_x)
        E_bot[j] += magz[Int(4*j)]
    end;



    res=string("flipp_up", "  ", "flipp_bot","  ", "mag","  ", "E_up","  ", "E_bot","\n")

    open("result.cvs", "a") do io
        write(io, res)
    end


    open("result.cvs", "a") do io
        writedlm(io, [flipp_up flipp_bot mag E_up E_bot])
    end


end


dostuff(parse(Int, ARGS[1]) , parse(Float64, ARGS[2]))
