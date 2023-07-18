using LinearAlgebra, SparseArrays, StaticArrays 
using StartUpDG
using Plots
using OrdinaryDiffEq
using SummationByPartsOperators
using MAT
using Trixi: AliveCallback

function initial_condition(x)

    if (0.4 <= x <= 0.6)
        b = (0.25*(1.0 + cos(10*pi*(x-0.5))))
    else
        b = 0.0
    end
    
    xl = 0.25
    xr = 0.75
    sigma_0 = 0.2

    if (xl <= x <= xr)
        A = (1 - sigma_0*(1 + cos(2*pi*(x - (xl + xr)/2.0)/(xr - xl)))) 
    else
        A=1.0
    end 


    h = 1.0 - b    
    u = 0.0

    return SVector(h * A, h * A * u, A, b)
end

function flux_ec(u_l, u_r)
    g = 1.0
    hA_L, hAu_L, A_L, b_L = u_l
    hA_R, hAu_R, A_R, b_R = u_r
    u_L = hAu_L / hA_L
    u_R = hAu_R / hA_R
    h_L, h_R = hA_L / A_L, hA_R / A_R

    hAu_avg = 0.5 * (hAu_L + hAu_R)
    u_avg = 0.5 * (u_L + u_R)

    hA = hAu_avg
    hAu = hAu_avg * u_avg + 0.5 * g * (A_L * h_L * h_R) + 0.5 * g * (A_L * h_L * b_R) 

    return SVector(hA, hAu, 0, 0)
end

# non-well-balanced version
function LxF_penalty(u_l, u_r)
    g = 1.0
    hA_L, hAu_L, A_L, b_L = u_l
    hA_R, hAu_R, A_R, b_R = u_r
    u_L = hAu_L / hA_L
    u_R = hAu_R / hA_R
    h_L, h_R = hA_L / A_L, hA_R / A_R
    
    c_L = abs(u_L) + sqrt(g * h_L)
    c_R = abs(u_R) + sqrt(g * h_R)
    lambda = max(c_L, c_R) 

    hA = 0.5 * lambda * (hA_R - hA_L)
    hAu = 0.5 * lambda * (hAu_R - hAu_L)

    return SVector(hA, hAu, 0, 0)
end


# well-balanced version with jumps of entropy variables
function LxF_penalty_wb(u_l, u_r)
    g = 1.0

    hA_L, hAu_L, A_L, b_L = u_l
    hA_R, hAu_R, A_R, b_R = u_r
    u_L = hAu_L / hA_L
    u_R = hAu_R / hA_R
    h_L, h_R = hA_L / A_L, hA_R / A_R
    c_L = abs(u_L) + sqrt(g * h_L)
    c_R = abs(u_R) + sqrt(g * h_R)
    lambda = max(c_L, c_R) 

    u = 0.5 * (u_L + u_R)
    h = 0.5 * (h_L + h_R)
    hA = 0.5 * (hA_L + hA_R)
    R = SMatrix{4, 4}([g * h + u^2 -u 0 0; -u 1 0 0; 0 0 0 0; 0 0 0 0] / hA)
    v_l = consvars2entropy(u_l)
    v_r = consvars2entropy(u_r)
    return 0.5 * lambda * R * (v_r - v_l)

end

function psi(u)
    g = 1.0
    hA, hAu, A, b = u
    h = hA / A
    u = hAu / hA
    return 0.5 * g * A * h * (h+b) * u
end


function rhs!(du::Matrix{<:SVector}, u, parameters, t)
    (; rd, md, D_skew) = parameters
    (; nxJ, J) = md

    fill!(du, zero(eltype(du)))

    uf = u[rd.Fmask, :]
    uP = uf[md.mapP]
    # for i in md.mapB
    #     hA, hAu, A, b = uf[i]
    #     uP[i] = SVector(hA, -hAu, A, b)
    # end
    interface_flux = @. flux_ec(uf, uP) * nxJ - LxF_penalty_wb(uf, uP)
    mul!(du, rd.LIFT, interface_flux)
    for e in 1:md.num_elements
        for i in axes(u, 1), j in axes(u, 1)
            du[i, e] += D_skew[i, j] * flux_ec(u[i, e], u[j, e])
        end
    end
    du ./= -J
end

function consvars2entropy(u)
    g = 1.0
    hA, hAu, A, b = u
    h, u = hA / A, hAu / hA
    return SVector(g * (h+b) - 0.5 * u^2, u, 0, 0)
end

N = 4
K = 200
rd = RefElemData(Line(), SBP(), N)

(VX,), EToV = uniform_mesh(Line(), K)
@. VX = 0.5 * (1 + VX) 
md = MeshData(((VX,), EToV), rd)


#md = MeshData(uniform_mesh(Line(), 2000), rd)
md = make_periodic(md)

Qr = rd.M * rd.Dr
D_skew = rd.M \ (Qr - Qr')
params = (; rd, md, D_skew)

rq, wq = gauss_quad(0, 0, N)
Vq = vandermonde(Line(), N, rq) / rd.VDM
Pq = (Vq' * diagm(wq) * Vq) \ (Vq' * diagm(wq))
xq = Vq * md.x
u = Pq * initial_condition.(xq)

du = similar(u)
rhs!(du, u, params, 0.0)
@show sum(dot.(consvars2entropy.(u), md.wJq .* du))

tspan = (0, 1.0)
ode = ODEProblem(rhs!, u, tspan, params)
h = estimate_h(rd, md)
sol = solve(ode, RK4(), adaptive=false, dt = 0.01 * h, callback=AliveCallback(alive_interval=10), 
saveat=LinRange(tspan..., 50), abstol=1e-12, reltol=1e-12)
u = sol.u[end]
u_0 = getindex.(sol.u[1], 1)
A_0 = getindex.(sol.u[1], 3)
h_0 = u_0./A_0
hAu_0 = getindex.(sol.u[1], 2)
hu_0 = hAu_0./A_0

hA = getindex.(u, 1)
hAu = getindex.(u, 2)
A = getindex.(u, 3)
b = getindex.(u, 4)

sol_md   = "wb_test_" * string(N)*"_K"*string(K)*".mat"
matwrite(sol_md, Dict( "h_0"=> h_0, "hu_0"=> hu_0, "Vp"=> rd.Vp, "wJq" => md.wJq, "x" => md.x, "hA" => hA, "hAu" => hAu, "A" => A, "b"=>b); compress = true)


