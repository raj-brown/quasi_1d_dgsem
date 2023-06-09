using MAT
using Plots
using LaTeXStrings
plot_font = "Computer Modern"
default(fontfamily=plot_font, framestyle=:box, label=nothing, grid=false)
scalefontsizes(1.0)

function initial_condition_rho(x)
    if x <= 0
        rho = 3.4718
    else
        rho  = 2.0
    end
    return rho
end


function initial_condition_u(x)
    
    if x < 0
     u = -2.5923 
    else
        u = -3.0
    end
    
    return u
end


function initial_condition_p(x)
    
    if x < 0
      p = 5.7118
    else
        p= 2.639
    end
    
    return p
end

function initial_condition(x)
    
    if x < 0
        rho, u, p, A = 3.4718, -2.5923, 5.7118, 1
    else
        rho, u, p, A = 2, -3, 2.639, 1.5
    end
    q = SVector(rho, u, p)
    return SVector(A * prim2cons(q, equations)..., A)
end



# Entropy Data
ec_data = matread("ec_t.mat")
N_skip =100

### Extract analytical data
t = vec(ec_data["t"])
ec = vec(ec_data["ec_val"])
x = ec_data["x"]
x = [x[1, :]; x[end, end]]

t = t[1:N_skip:end ]
ec = ec[1:N_skip:end ]


##Initial values
rho = initial_condition_rho.(x)
u = initial_condition_u.(x)
p = initial_condition_p.(x)



#plot(t, ec, linewidth=3, c=:red, legend=:false)
#ylims!((10^-15, 10^0))
#plot!( yscale=:log10, minorgrid=true)

#xlabel!("Time (t)")
#ylabel!("Entropy residuals")
#savefig("ec.png") 

plot(x, rho, linewidth=1.5, c=:black, legend=:bottomleft, labels="ρ" )
plot!(x, u, linewidth=1.5, c=:red, legend=:bottomleft, labels="u" )
plot!(x, p, linewidth=1.5, c=:blue, legend=:bottomleft, labels="p" )
ylims!((-5, 6))

xlabel!("x ")
savefig("rho_ic.png") 
