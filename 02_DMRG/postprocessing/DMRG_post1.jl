using DataFrames
using Plots
using LaTeXStrings
using CSV

## READ Data
file="/home/phil/develop/SquareIce/06_data/01_DMRG/20201213_DMRG_ED_comparison_EA_EB_Entropy/dmrg.dat"
df=CSV.File(file)


## %% Plot chess operators
p=[plot(),plot(),plot(),plot()]
gr() # Set the backend to Plotly

for L = unique(df.Number_Plaquettes)
    println("plotting chess L=", L)
    Index=findall(df.Number_Plaquettes.==L)
    Oflip=[(df[i].coupling,df[i].Oflip) for i in Index]
    Oflipp=[(df[i].coupling,df[i].Oflipp) for i in Index]
    winding=[(df[i].coupling,df[i].winding_number) for i in Index]
    E=[(df[i].coupling,df[i].Energy_GS) for i in Index]

    plot!(p[1],Oflip, label=string("L=",string(L)),marker="+")
    plot!(p[2],Oflipp, label=string("L=",string(L)), marker="x")
    plot!(p[3],winding, label=string("L=",string(L)), marker="x")
    plot!(p[4],E, label=string("L=",string(L)), marker="x")
end
xlabel!("\$\\lambda\$")
ylabel!(p[1],"\$<{O}_\\mathrm{flip}>\$")
ylabel!(p[2],"\$<{O}_\\mathrm{flipp}>\$")
ylabel!(p[3],"\$<{W}_y>\$")
ylabel!(p[3],"\$<E>\$")
#ylabel!(p[3],"\$<{M}_A^2+{M}_B^2>\$")
plot(p[1],p[2],p[3],p[4],layout=(4,1), legend=:bottomleft)
#current()
# save figure
savefig("01_notes/images/20201214_Oflip.pdf")

p=plot()
gr()
for k =1:10
    entropy = (0:df[k].Number_Plaquettes,[parse(Float64,x) for x in split(df[k].entropy[2:end-1],",")])
    plot!(p,entropy,label=string("lambda=",string(df[k].coupling)))
end
xlabel!("\$x\$")
ylabel!("entropy ")
plot(p, legend=:bottomleft)
current()
savefig("01_notes/images/20201214_entropy.pdf")
