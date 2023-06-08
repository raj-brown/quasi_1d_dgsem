#=
This code test the performs the convergence study for 1d-ESDG method for with refrence solution constructed 
    by using manufactured solution methos
=#


using LinearAlgebra, SparseArrays, StaticArrays 
using StartUpDG
using Plots
using Revise
using OrdinaryDiffEq
using Trixi: cons2prim, prim2cons, cons2entropy, ln_mean, inv_ln_mean, CompressibleEulerEquations1D, DissipationLocalLaxFriedrichs
using MAT
using ForwardDiff

function primitive_manufactured_sol()
    A(x, t) =   1.0 - 0.1 * (1.0 + cos(pi * (x - 0.5) / 0.5))
    rho(x, t) = (1.0 + 0.1*sin(2*pi*x) + 0.1*cos(2*pi*x))*exp(-t)
    u(x, t) =   (1.0 + 0.1*sin(2*pi*x) + 0.1*cos(2*pi*x))*exp(-t)
    p(x, t) =   (1.0 + 0.1*sin(2*pi*x) + 0.1*cos(2*pi*x))*exp(-t)    
    return A, rho, u, p
end

function source(x, t, equations)
    A, rho, u, p = primitive_manufactured_sol()
    E(x, t) = (1.0 / 2.0)*rho(x, t) * u(x, t)^2 + p(x,t)/(equations.gamma - 1.0)
    
    f1(x, t) = ForwardDiff.derivative(t -> A(x,t)*rho(x,t), t) + ForwardDiff.derivative(x -> A(x,t) * rho(x,t) * u(x,t), x) 
    
    f2(x, t) = ForwardDiff.derivative(t -> A(x,t)*rho(x,t) * u(x, t), t) + ForwardDiff.derivative(x -> A(x, t)*rho(x, t) *u(x, t) * u(x, t), x) + A(x, t) * ForwardDiff.derivative(x -> p(x, t), x)
    
    f3(x, t) = ForwardDiff.derivative(t -> A(x, t)*E(x, t), t) + ForwardDiff.derivative( x -> A(x, t)*u(x, t)*(E(x, t) + p(x, t)), x) 

    return SVector(f1(x, t), f2(x, t), f3(x, t), 0.0)
end

function manufacture_sol(x0, t0, equations::CompressibleEulerEquations1D)
    # generate functions
    A, rho, u, p = primitive_manufactured_sol()
    rho_0 = rho(x0, t0)
    u_0 =  u(x0, t0)
    p_0 =  p(x0, t0)
    q = SVector(rho_0, u_0, p_0)
    return SVector(A(x0, t0) * prim2cons(q, equations)..., A(x0, t0)) 
end 

function A2cons(uA, ::CompressibleEulerEquations1D)
    A = uA[4]
    return SVector{3}(uA[1:3] ./ A)
end

function flux_ec(uA_ll, uA_rr, equations::CompressibleEulerEquations1D)
    gamma = equations.gamma
    u_ll = A2cons(uA_ll, equations)
    u_rr = A2cons(uA_rr, equations)
    A_ll, A_rr = uA_ll[4], uA_rr[4]
    
    rho_ll, v1_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, p_rr = cons2prim(u_rr, equations)
    
    rho_mean = ln_mean(rho_ll, rho_rr)
    inv_rho_p_mean = p_ll * p_rr * inv_ln_mean(rho_ll * p_rr, rho_rr * p_ll)
    v1_avg = 0.5 * (v1_ll + v1_rr)
    p_avg  = 0.5 * (p_ll + p_rr)
    velocity_square_avg = 0.5 * (v1_ll * v1_rr)
    
    v1A_ll = v1_ll * A_ll 
    v1A_rr = v1_rr * A_rr
    v1A_avg = 0.5 * (v1A_ll + v1A_rr)
    
    f1 = rho_mean * v1A_avg
    f2 = f1 * v1_avg + A_ll * p_avg
    f3 = f1 * ( velocity_square_avg + inv_rho_p_mean / (gamma-1) ) + 0.5 * (p_ll * v1A_rr + p_rr * v1A_ll)
    
    return SVector(f1, f2, f3, 0.0)
end

LxF_dissipation = DissipationLocalLaxFriedrichs()

function LxF_penalty(u_l, u_r, equations::CompressibleEulerEquations1D)
    val = -LxF_dissipation(u_l, u_r, 1, equations) # this is -lambda / 2 * (u_r - u_l)
    return SVector(val[1], val[2], val[3], 0)
end

function rhs!(du::Matrix{<:SVector}, u, parameters, t)
    (; rd, md, D_skew, equations, source) = parameters
    (; nxJ, J) = md

    fill!(du, zero(eltype(du)))

    uf = u[rd.Fmask, :]
    uP = uf[md.mapP]
    # for i in md.mapB
    #     rhoA, rhoAu, E, A = uf[i]
    #     uP[i] = SVector(rhoA, -rhoAu, E, A)
    #     # uP[i] = initial_condition(md.xf[i], equations)
    # end
    interface_flux = @. flux_ec(uf, uP, equations) * nxJ - LxF_penalty(uf, uP, equations)
    mul!(du, rd.LIFT, interface_flux)
    for e in 1:md.num_elements
        for i in axes(u, 1), j in axes(u, 1)
            du[i, e] += D_skew[i, j] * flux_ec(u[i, e], u[j, e], equations)
        end
    end
    du ./= -J 
    @. du = du + source(md.x, t, equations)
end

println("Started Convergence!")  

N_arr = [1, 2, 3, 4, 5]
K_arr = [2, 2^2, 2^3, 2^4, 2^5]

for i =eachindex(N_arr)
    for j =eachindex(N_arr)
        println("Started Convergence!")     
        N = N_arr[i]
        K = K_arr[j]
        rd = RefElemData(Line(), SBP(), N)
        md = MeshData(uniform_mesh(Line(), K), rd)
        md = make_periodic(md)
        equations = CompressibleEulerEquations1D(1.4)
        Qr = rd.M * rd.Dr
        D_skew = rd.M \ (Qr - Qr')
        params = (; rd, md, D_skew, equations, source)
        println("Initializing Manufacture Solution")
        u = manufacture_sol.(md.x, 0.0, equations) 

        # check entropy residual 
        println("Checking Entropy Residuals!")
        du = similar(u)
        rhs!(du, u, params, 0.0)
        w = map(x -> SVector{4}(x..., 0.0), cons2entropy.(A2cons.(u, equations), equations))
        @show sum(dot.(w, md.wJq .* du))

        tspan = (0, .1)
        ode = ODEProblem(rhs!, u, tspan, params)
        println("Computing...")
        sol = solve(ode, RK4(), saveat=LinRange(tspan..., 50), abstol=1e-10, reltol=1e-10)
        sol_md = "Num_Manu_sol_N"*string(N)*"_K"*string(K)*".mat"
        sol_md_a = "Ana_Manu_sol_N"*string(N)*"_K"*string(K)*".mat"

        u = sol.u[end]
        rhoA = getindex.(u, 1)
        rhouA = getindex.(u, 2)
        EA = getindex.(u, 3)
        A = getindex.(u, 4)
 

        u_a = manufacture_sol.(md.x, tspan[end], equations)
        rhoA_a = getindex.(u_a, 1)
        rhouA_a = getindex.(u_a, 2)
        EA_a = getindex.(u_a, 3)
        A_a = getindex.(u_a, 4)
       
          
        matwrite(sol_md, Dict( "wJq" => md.wJq, "x" => md.x, "rho" => rhoA ./ A, "rhou" => rhouA ./ A, "E"=>EA ./ A); compress = true)
        matwrite(sol_md_a, Dict( "wJq" => md.wJq, "x" => md.x, "rho_a" => rhoA_a ./ A_a, "rhou_a" => rhouA_a ./ A_a, 
        "E_a"=>EA_a ./ A_a); compress = true)
    end
end

