
#this is a set of funtions to generate every possible combinations of a set of N element. The values that our vector can have are stored in an array (in our case [-1,1])
#I use [-1,1] out of simplicisty to compure the gauss law expectatio value
import Base.iterate,Base.length
struct Combinations{T}
    itr::Vector{T}
    count::Int64
    itrsize::Int64
    function Combinations(itr::Vector{T},count::Int) where T
        new{T}(itr,Int64(count),length(itr))
    end
end

function iterate(c::Combinations,state::Int64=0)
    if state>=length(c)
        return nothing
    end
    indices=digits(state,base=c.itrsize,pad=c.count)
    [c.itr[i] for i in (indices .+1)],state+1
end

function length(c::Combinations)
    length(c.itr) ^ c.count
end

function iterate(c::Combinations,state::Int64=0)
    if state>=length(c)
        return nothing
    end
    indices=digits(state,base=c.itrsize,pad=c.count)
    [c.itr[i] for i in (indices .+1)],state+1
end

# this function maps a vector of lenght (2*x-1)*y) into a rectangular lattice 
#I have 2*x-1 spin in x dirextion (I call x what we usually call Lx ) because i treat the orizontal and vertical spin equally. 


function rectangularize(x,y)
    linear_combination=collect(Combinations([-1,1],(2*x-1)*y))
    rectangulare_lattice_cell=zeros(2*x-1,y)
    rectangular_lattice=zeros(2^((2*x-1)*y),2*x-1,y)

    for k in [1:1:2^((2*x-1)*y);];
        for j in [1:1:y;];
            for i in [1:1:2*x-1;];
                rectangular_lattice[k,i,j]=linear_combination[k][i+j-1]
            end
        end
    end
    return rectangular_lattice
end

#given a rectangular lattice, returns the expectation value of the gauss operator starting from the vertical spin 
#(work onyly in the odd sites for my setup) going into the negative direction into the y axis tG
function gauge_invariant_site_exp_val(array,x,y,x_max,y_max)
    expectation_value_G=0
    if  x< 2*x_max -1  && x> 1 && y> 1 
        expectation_value_G=array[x,y]+array[x,y-1]+array[x-1,y]+array[x+1,y] 
    elseif x< 2*x_max -1 && x> 1 && y==1
        expectation_value_G=array[x,y]+array[x,y_max]+array[x-1,y]+array[x+1,y] 
    elseif x== 1   && y> 1    
        expectation_value_G=array[x,y]+array[x,y-1]+(-1)^(y+1)+array[x+1,y] 
    elseif x== 1   && y==1    
        expectation_value_G=array[x,y]+array[x,y_max]+(-1)^(y+1)+array[x+1,y] 
    elseif x== 2*x_max -1  && y> 1   
        expectation_value_G=array[x,y]+array[x,y-1]+array[x-1,y]+(-1)^(y+1) 
    elseif x== 2*x_max -1  && y==1
        expectation_value_G=array[x,y]+array[x,y_max]+array[x-1,y]+(-1)^(y+1) 
    end
    return expectation_value_G    
end


#this function loops over all possible states (enumerated with index k) generate in rectangularize and for evey single state check that the expectation value of the Gauus operator
#in every nod is zero. It True, add the index k the the array states and return states

function counting_gauge_invariant_states(array,x_max,y_max)
    check=true
    states=[]
    
    for k in [1:1:size(array[:,1,1])[1];];
        i=1
        j=1
        while i < 2*x_max && check;
            while j < y_max+1 && check; 
                check =  abs(gauge_invariant_site_exp_val(array[k,:,:] ,i,j,x_max,y_max)) > 0.0001
                j+=1
            end
            j=1
            i+=1
        end
        if check push!(states, k) end
        check=true
    end

return states
end

Lx=5  
Ly=2

array=rectangularize(Lx,Ly)
states=counting_gauge_invariant_states(array,Lx,Ly)
num_states=size(states)[1]
println("Lx=$Lx   Ly=$Ly   Number of states = $num_states " )