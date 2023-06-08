using LinearAlgebra, SparseArrays
using StaticArrays 
using StartUpDG
using Plots
using OrdinaryDiffEq
using Trixi: AliveCallback, pressure, cons2prim, prim2cons, cons2entropy, ln_mean, inv_ln_mean, CompressibleEulerEquations1D, DissipationLocalLaxFriedrichs
using MAT
function a(x)
    x = abs(x)
    if (0.0 <= x <= 5.0)
        return 1 + 1.5 * (1 - x / 5.0)^2
    elseif (5.0 <= x <= 10.0)
        return 1 + 0.5 * (1 - x / 5.0)^2
    else
        println("$x value out of range!")
    end
end

p_inflow() = 9.608491914104024e04
 
p_outflow() = 8.497381834742936e04

@show typeof(p_inflow)

function initial_condition(x, equations::CompressibleEulerEquations1D) 

    T0 = 300.0 # Kelvin
    R = 287.0
    Ma = 0.2395 # Mach number = u / c

    rho_inflow = p_inflow() / (R * T0)
    c = sqrt(equations.gamma * p_inflow() / rho_inflow)
    u_inflow = Ma * c

    A = a(x)
    q = SVector(rho_inflow, u_inflow, p_inflow())
    return SVector(A * prim2cons(q, equations)..., A)
end


function A2cons(uA, ::CompressibleEulerEquations1D)
    A = uA[4]
    return SVector{3}(uA[1:3]./A)
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

    # inflow: set rho and p (two thermodynamic variables)
    rhoA, rhoAu, _, A = uP[1]
    rho_inflow, u_inflow, _ = cons2prim(A2cons(uP[1], equations), equations)    
    rhoA, rhoAu, _ = initial_condition(0, equations)
    rho, rho_v1, rho_e = prim2cons(SVector(rho_inflow, u_inflow, p_inflow()), equations)
    uP[1] = SVector(rhoA, rhoA * u_inflow, A * rho_e, A) # left BC

    # outflow    
    rhoA, rhoAu, _, A = uP[end]
    rho_outflow, u_outflow, _ = cons2prim(A2cons(uP[end], equations), equations)    
    rho, rho_v1, rho_e = prim2cons(SVector(rho_outflow, u_outflow, p_outflow()), equations)
    uP[end] = SVector(A * rho, A * rho_v1, A * rho_e, A) # Right BC

    interface_flux = @. flux_ec(uf, uP, equations) * nxJ - LxF_penalty(uf, uP, equations)
    mul!(du, rd.LIFT, interface_flux)
    for e in 1:md.num_elements
        for i in axes(u, 1), j in axes(u, 1)
            du[i, e] += D_skew[i, j] * flux_ec(u[i, e], u[j, e], equations)
        end
    end
    du ./= -J
end

N = 3
K = 64

rd = RefElemData(Line(), SBP(), N)
(VX,), EToV = uniform_mesh(Line(), K)
@. VX = 0.5 * (1 + VX) * 10

md = MeshData(((VX,), EToV), rd)
equations = CompressibleEulerEquations1D(1.4)
Qr = rd.M * rd.Dr
D_skew = rd.M \ (Qr - Qr')
params = (; rd, md, D_skew, equations)
u = initial_condition.(md.x, equations)
     
println("Checking Entropy Residuals!")
du = similar(u)
rhs!(du, u, params, .1)
w = map(x -> SVector{4}(x..., 0.0), cons2entropy.(A2cons.(u, equations), equations))
@show sum(dot.(w, md.wJq .* du))

println("Computing the ODE solution...")
tspan = (0, 1.0)
ode = ODEProblem(rhs!, u, tspan, params)
sol = solve(ode,  RDPK3SpFSAL49(), saveat=LinRange(tspan..., 100), 
            abstol=1e-10, reltol=1e-10, callback=AliveCallback(alive_interval=100))

@gif for i in eachindex(sol.u)
    u = sol.u[i]
    rhoA = getindex.(u, 1)
    rhouA = getindex.(u, 2)
    A = getindex.(u, 4)

    rho = rhoA ./ A
    p = pressure.(A2cons.(u, equations), equations)
    v = rhouA ./ rhoA
    c = @. sqrt(equations.gamma * p / rho)
    Ma = @. rhouA / (rhoA * c)
    #plot(vec(rd.Vp * md.x), vec(rd.Vp * Ma), leg=false, ylim = (0, .8))
    plot(vec(rd.Vp * md.x), vec(rd.Vp * p), leg=false)
    title!("Time t=$(sol.t[i])")
end

u = sol.u[end]
p = pressure.(A2cons.(u, equations), equations)
x_p = rd.Vp * md.x
p_p = rd.Vp * p



u = sol.u[end]
rhoA = getindex.(u, 1)
rhouA = getindex.(u, 2)
A = getindex.(u, 4)

rho = rhoA ./ A
p = pressure.(A2cons.(u, equations), equations)
c = @. sqrt(equations.gamma * p / rho)
Ma = @. rhouA / (rhoA .* c)
Ma = rd.Vp * Ma
plot(rd.Vp * md.x, rd.Vp * p, leg=false)
title!("Time t=$(sol.t[end])")

matwrite("trans_subsonic.mat", Dict( "x" => x_p, "pressure" => p_p, "Ma" => Ma); compress = true)




