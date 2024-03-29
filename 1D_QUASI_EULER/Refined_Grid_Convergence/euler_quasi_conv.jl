using LinearAlgebra, SparseArrays, StaticArrays 
using StartUpDG
using Plots
using Revise
using OrdinaryDiffEq
using Trixi: cons2prim, prim2cons, cons2entropy, ln_mean, inv_ln_mean, CompressibleEulerEquations1D, DissipationLocalLaxFriedrichs
using MAT

using ForwardDiff

function initial_condition(x, equations::CompressibleEulerEquations1D) 
    A = 1.0 - 0.2 * (1 + cos(pi * (x - 0.5) / 0.5))
    rho = 1.0 - 0.1 * (1 + sin(pi * (x - 0.1) / 0.5))
    u   = 0.0
    p   = rho^equations.gamma
    q = SVector(rho, u, p)
    return SVector(A * prim2cons(q, equations)..., A)
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
    (; rd, md, D_skew, equations) = parameters
    (; nxJ, J) = md

    fill!(du, zero(eltype(du)))

    uf = u[rd.Fmask, :]
    uP = uf[md.mapP]
    interface_flux = @. flux_ec(uf, uP, equations) * nxJ - LxF_penalty(uf, uP, equations)
    mul!(du, rd.LIFT, interface_flux)
    for e in 1:md.num_elements
        for i in axes(u, 1), j in axes(u, 1)
            du[i, e] += D_skew[i, j] * flux_ec(u[i, e], u[j, e], equations)
        end
    end
    du ./= -J
end

N_arr = [1, 2, 3, 4, 5]
K_arr = [1, 2, 3, 4, 5]

for i = eachindex(N_arr)
    for j =eachindex(N_arr)
        println("Started Convergence!")     
        N_loc = N_arr[i]
        K_loc = 2^K_arr[j]
        rd = RefElemData(Line(), SBP(), N_loc)
        md = MeshData(uniform_mesh(Line(), K_loc), rd)
        md = make_periodic(md)
        equations = CompressibleEulerEquations1D(1.4)
        Qr = rd.M * rd.Dr
        D_skew = rd.M \ (Qr - Qr')
        params = (; rd, md, D_skew, equations)
        u = initial_condition.(md.x, equations)     
        # check entropy residual 
        println("Checking Entropy Residuals!")
        du = similar(u)
        rhs!(du, u, params, 0.0)
        w = map(x -> SVector{4}(x..., 0.0), cons2entropy.(A2cons.(u, equations), equations))
        @show sum(dot.(w, md.wJq .* du))
        tspan = (0, 0.1)
        ode = ODEProblem(rhs!, u, tspan, params)
        println("Computing...")
        sol = solve(ode, RK4(), saveat=LinRange(tspan..., 50), abstol=1e-10, reltol=1e-10)
        sol_md = "sol_N"*string(N_loc)*"_K"*string(K_loc)*".mat"
        u = sol.u[end]
        rhoA = getindex.(u, 1)
        rhouA = getindex.(u, 2)
        EA = getindex.(u, 3)

        A = getindex.(u, 4)
        matwrite(sol_md, Dict( "wJq" => md.wJq, "x" => md.x, "rho" => rhoA ./ A, "rhou" => rhouA ./ A, "E"=>EA ./ A); compress = true)
        
    end
end

