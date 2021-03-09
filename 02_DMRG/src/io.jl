using JLD2
using FileIO

function save_configuration( A, Lx, BondD, lambda; theta = 0.0 , mu = 0.0, Ly = 2, prefix = "SQ")

    fname = string(prefix , "_Lx-", Lx, "_Ly-", Ly, "_BondD-",BondD, "_lambda-",lambda,
            "_theta-",theta, "_mux-", mu[1], "_muy-", mu[2], ".jld2" )
    print("Loading: ", fname, "\n")
    print("Configuration saved in: ", fname ,"\n")
    save(fname, "MPS", A , "Lx", Lx, "D" , BondD , "lambda" , lambda , "theta" , theta , "mu", mu, "Ly", Ly)

end

function load_configuration( Lx, BondD, lambda; theta = 0 , mu = 0, Ly = 2, prefix = "SQ")

    fname = string(prefix , "_Lx-", Lx, "_Ly-", Ly, "_BondD-",BondD, "_lambda-",lambda,
            "_theta-",theta, "_mux-", mu[1], "_muy-", mu[2], ".jld2" )
    print("Loading: ", fname, "\n")
    dict = load(fname)

    return dict
end

function load_configuration( fname)

    print("Loading: ", fname, "\n")
    dict = load(fname)
    return dict
end
