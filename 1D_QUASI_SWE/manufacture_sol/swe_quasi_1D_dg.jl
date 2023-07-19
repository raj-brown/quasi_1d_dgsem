using LinearAlgebra, SparseArrays, StaticArrays 
using StartUpDG
using Plots
using OrdinaryDiffEq
using SummationByPartsOperators
using MAT
using Trixi: AliveCallback
using ForwardDiff

function primitive_manufactured_solution()   
    A(x, t) = exp(sin(2*pi*x))
    b(x, t) = sin(pi*x) * sin(pi*x)
    h(x, t) = 3 + 0.1*exp(cos(2*pi*x))*exp(-t)
    u(x, t) = (sin(cos(2*pi*x)))*exp(-t)
    return A, b, h, u
end

function source(x, t)
    A, b, h, u = primitive_manufactured_solution()
    g = 1.0
    f1(x, t) = ForwardDiff.derivative(t-> A(x, t)*h(x, t), t) + 
               ForwardDiff.derivative(x->A(x, t)*h(x, t)*u(x, t), x)
    f2(x, t) = ForwardDiff.derivative(t-> A(x, t)*h(x, t)*u(x, t), t) + 
               ForwardDiff.derivative(x->A(x,t)*h(x, t) * u(x, t)^2, x) + 
               g*A(x,t)*h(x, t) * ForwardDiff.derivative(x->h(x, t) + b(x, t), x)

    return SVector(f1(x, t), f2(x, t), 0.0, 0.0)
   
end

function initial_condition(x)  
    A, b, h, u = primitive_manufactured_solution()
    A_0 = A(x, 0)
    b_0 = b(x, 0)
    h_0 = h(x, 0)
    u_0 = u(x, 0)
    return SVector(h_0* A_0, h_0 * A_0 * u_0, A_0, b_0)
end


function analytical_sol(x, t)  
    A, b, h, u = primitive_manufactured_solution()
    A_t = A(x, t)
    b_t = b(x, t)
    h_t = h(x, t)
    u_t = u(x, t)
    return SVector(h_t* A_t, h_t * A_t * u_t, A_t, b_t)
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
    
    interface_flux = @. flux_ec(uf, uP) * nxJ - LxF_penalty(uf, uP)
    mul!(du, rd.LIFT, interface_flux)
    for e in 1:md.num_elements
        for i in axes(u, 1), j in axes(u, 1)
            du[i, e] += D_skew[i, j] * flux_ec(u[i, e], u[j, e])
        end
    end
    du ./= -J
    @. du = du + source(md.x, t)
end

function consvars2entropy(u)
    g = 1.0
    hA, hAu, A, b = u
    h, u = hA / A, hAu / hA
    return SVector(g * (h+b) - 0.5 * u^2, u, 0, 0)
end

N_arr = [1, 2, 3, 4, 5]
K_arr = [16, 32, 64, 128, 256]
#K_arr = [2, 4, 8, 16, 32]


for i =eachindex(N_arr)
    for j =eachindex(N_arr)
        println("Started Convergence!")     
        N = N_arr[i]
        K = K_arr[j]

        rd = RefElemData(Line(), SBP(), N)

        (VX,), EToV = uniform_mesh(Line(), K)
        #@. VX = 0.5 * (1 + VX) 
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

        tspan = (0, 0.1)
        ode = ODEProblem(rhs!, u, tspan, params)
        h = estimate_h(rd, md)
        sol = solve(ode, RK4(), adaptive=true, dt = 0.01 * h, callback=AliveCallback(alive_interval=10), 
        saveat=LinRange(tspan..., 50), abstol=1e-12, reltol=1e-12)


        sol_md = "Num_Manu_sol_N"*string(N)*"_K"*string(K)*".mat"
        sol_md_a = "Ana_Manu_sol_N"*string(N)*"_K"*string(K)*".mat"
        u = sol.u[end]
        hA = getindex.(u, 1)
        hAu = getindex.(u, 2)
        A = getindex.(u, 3)
        b = getindex.(u, 4)

        u_a = analytical_sol.(md.x, tspan[end])
        hA_a = getindex.(u_a, 1)
        hAu_a = getindex.(u_a, 2)
        A_a = getindex.(u_a, 3)
        b_a = getindex.(u_a, 4)

        matwrite(sol_md, Dict( "wJq" => md.wJq, "x" => md.x, "hA" => hA, "hAu" => hAu, "A" => A, 
        "b"=>b); compress = true)
        matwrite(sol_md_a, Dict( "wJq" => md.wJq, "x" => md.x, "hA_a" => hA_a, "hAu_a" => hAu_a, "A_a" => A_a,
        "b_a"=>b_a); compress = true)

        
    end
end


# @gif for u in sol.u
#     hA = getindex.(u, 1)
#     A = getindex.(u, 3)
#     b = getindex.(u, 4)
#     plot(rd.Vp * md.x, rd.Vp * (b + hA ./ A))
#     plot!(rd.Vp * md.x, rd.Vp * b, leg=false)
# end