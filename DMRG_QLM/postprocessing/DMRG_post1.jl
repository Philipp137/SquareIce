using DataFrames
using Plots
using LaTeXStrings

## READ Data
dir="/home/phil/develop/SquareIce/big_fss/"
#data=readtable(fname, separator=' ')
data=[]
lambda=0;Lsize=0
counter = 0
for (root, dirs, files) in walkdir(dir)
    global data
    global df
    global counter
    for file in files
        fpath=joinpath.(root, file)
        if file == "scaling" && filesize(fpath)>0 # read only nonempty file called scaling
            println("reading: ", fpath)
            ## read data from file
            df=readtable(fpath, separator=' ')
            ## insert lattice size
            m=match(r"/L(?<Lsize>\d+)",root)
            insert!(df,1,parse(Int64,m[:Lsize]),:Lsize)
            ## append the data or create new array
            if counter==0
                data =DataFrame(df[size(df,1),:])
            else
                data=push!(data,df[size(df,1),:])
                # data = last(df)
            end
            counter += 1
        else
            #println(df)
        end
    end
end
#return df


## %% Plot chess operators
p=[plot(),plot(),plot()]
gr() # Set the backend to Plotly
for L in unique(data.Lsize)
    println("plotting chess L=", L)
    Ma=data[data.Lsize.==L,:].coupling,data[data.Lsize.==L,:].chess_up
    Mb=data[data.Lsize.==L,:].coupling,data[data.Lsize.==L,:].chess_down
    plot!(p[1],Ma, label=string("L=",string(L)), ls=(:dot),marker = (:dot))
    plot!(p[2],Mb, label=string("L=",string(L)), ls=(:dot),marker = (:dot))
    if L==60
        plot!(p[3],Ma[1],Ma[2].^2 + Mb[2].^2,label=string("L=",string(L)), ls=(:dot),marker = (:dot,:hexagon))
    elseif L==20
        plot!(p[3],Ma[1],Ma[2].^2 + Mb[2].^2,label=string("L=",string(L)), ls=(:dashed),marker = (:hexagon))
    else

    end
end
xlabel!("\$\\lambda\$")
ylabel!(p[1],"\$<{M}_A>\$")
ylabel!(p[2],"\$<{M}_B>\$")
ylabel!(p[3],"\$<{M}_A^2+{M}_B^2>\$")
plot(p[1],p[2],p[3],layout=(3,1),legend=:bottomleft)
current()
# save figure
savefig(string( dir,"/../notes/images/chessop.pdf"))
