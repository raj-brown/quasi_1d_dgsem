using LinearAlgebra, SparseArrays, StaticArrays 
using StartUpDG
using Plots
using OrdinaryDiffEq
using SummationByPartsOperators
using MAT
using Trixi: AliveCallback

function initial_condition(x)
    xl = 100.0
    xr = 400.0
       
    if (xl <= x <= xr)
        A = 5.0 - 0.7065*(1.0 + cos(2*pi *( x - 250)/300))
    else
        A = 5.0
    end
    b = 0.0
    h = 2.0
    Q = 20.0
    H = A*h
    u = Q/H 
    return SVector(h * A, h * A * u, A, b)
end

function flux_ec(u_l, u_r)
    g = 9.81
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
    g = 9.81
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
    g = 9.81
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
    R = SMatrix{4, 4}([(g * h + u^2) -u 0 0; -u 1 0 0; 0 0 0 0; 0 0 0 0] / hA)
    v_l = consvars2entropy(u_l)
    v_r = consvars2entropy(u_r)
    return 0.5 * lambda * R * (v_r - v_l)

end

function psi(u)
    g = 9.81
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

    # Inflow boundary condition-
    hA, hAu, A, b = uP[1]
    h = hA/A    
    u_in = 20.0/(h*A)
    h_in = 2.0
    uP[1] = SVector(h_in*A, h_in*A*u_in, A, b)


    #  # Outflow boundary condition-
     hA, hAu, A, b = uP[end]
     h_out = 1.85
     u_out = hAu/hA
     uP[end] = SVector(h_out*A, h_out*A*u_out, A, b)

    interface_flux = @. flux_ec(uf, uP) * nxJ - LxF_penalty(uf, uP)
    mul!(du, rd.LIFT, interface_flux)
    for e in 1:md.num_elements
        for i in axes(u, 1), j in axes(u, 1)
            du[i, e] += D_skew[i, j] * flux_ec(u[i, e], u[j, e])
        end
    end
    du ./= -J
end

function consvars2entropy(u)
    g = 9.81
    hA, hAu, A, b = u
    h, u = hA / A, hAu / hA
    return SVector(g * (h+b) - 0.5 * u^2, u, 0, 0)
end

N = 2
K = 200
rd = RefElemData(Line(), SBP(), N)

(VX,), EToV = uniform_mesh(Line(), K)
@. VX = 0.5 * (1 + VX)*500 
md = MeshData(((VX,), EToV), rd)


#md = MeshData(uniform_mesh(Line(), 2000), rd)
#md = make_periodic(md)

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

tspan = (0, 10000)
ode = ODEProblem(rhs!, u, tspan, params)
h = estimate_h(rd, md)
sol = solve(ode, RK4(), adaptive=true, callback=AliveCallback(alive_interval=10),saveat=LinRange(tspan..., 50))
#sol = solve(ode, RDPK3SpFSAL49(), callback=AliveCallback(alive_interval=10), 
#saveat=LinRange(tspan..., 50), abstol=1e-10, reltol=1e-9)
u = sol.u[end]
hA = getindex.(u, 1)
hAu = getindex.(u, 2)
A = getindex.(u, 3)
b = getindex.(u, 4)
sol_md = "Sol_N" * string(N)*"_K"*string(K)*".mat"
matwrite(sol_md, Dict( "wJq" => md.wJq, "x" => md.x, "hA" => hA, "Vp"=> rd.Vp, "hAu" => hAu, "A" => A, "b"=>b); compress = true)
