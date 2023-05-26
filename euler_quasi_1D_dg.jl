using LinearAlgebra, SparseArrays, StaticArrays 
using StartUpDG
using Plots
using OrdinaryDiffEq
using Trixi: prim2cons, cons2entropy, ln_mean, inv_ln_mean, CompressibleEulerEquations1D

function s(x)
    if (x >= 0) & x <= 5
        f = 1 + 1.5*(1 - x/5)^2
    else
        if x >= 5 && x <= 10
            f = 1 + 0.5*(1 - x/5)^2;
        end 
    end
end


function initial_condition(x, equations::CompressibleEulerEquations1D)
    
    A = 1 - 0.2 * (1 + cos(pi * (x - 0.5) / 0.5))

    rho = 1.0 + .5 * exp(-100 * x^2)
    u   = 0.0
    p   = rho^equations.gamma

    # if x < 0
    #     rho, u, p, A = 3.4718, -2.5923, 5.7118, 1
    # else
    #     rho, u, p, A = 2, -3, 2.639, 1.5
    # end

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
    for i in md.mapB
        rhoA, rhoAu, E, A = uf[i]
        uP[i] = SVector(rhoA, -rhoAu, E, A)
        # uP[i] = initial_condition(md.xf[i], equations)
    end
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
rd = RefElemData(Line(), SBP(), N)
md = MeshData(uniform_mesh(Line(), 100), rd)
# md = make_periodic(md)

equations = CompressibleEulerEquations1D(1.4)

Qr = rd.M * rd.Dr
D_skew = rd.M \ (Qr - Qr')
params = (; rd, md, D_skew, equations)
u = initial_condition.(md.x, equations)

# rq, wq = gauss_quad(0, 0, N)
# Vq = vandermonde(Line(), N, rq) / rd.VDM
# Pq = (Vq' * diagm(wq) * Vq) \ (Vq' * diagm(wq))
# xq = Vq * md.x
# u = Pq * initial_condition.(xq, equations)

# check entropy residual 
du = similar(u)
rhs!(du, u, params, 0.0)
w = map(x -> SVector{4}(x..., 0.0), cons2entropy.(A2cons.(u, equations), equations))
@show sum(dot.(w, md.wJq .* du))

tspan = (0, 4.0)
ode = ODEProblem(rhs!, u, tspan, params)
sol = solve(ode, RK4(), saveat=LinRange(tspan..., 50))

@gif for u in sol.u
    rhoA = getindex.(u, 1)
    A = getindex.(u, 4)
    plot(rd.Vp * md.x, rd.Vp * (rhoA ./ A), leg=false, ylims=(.5, 1.5))
end

