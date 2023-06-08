using LinearAlgebra
using Plots
using MAT
using Revise

error_mat = zeros(5, 5)
N_arr = [1, 2, 3, 4, 5]
K_arr = [20, 40, 80, 160, 320]

for i=eachindex(N_arr)
    for j =eachindex(K_arr)
        N = i
        K = 2^j
        sol_md = "Num_Manu_sol_N"*string(N)*"_K"*string(K)*".mat"
        sol_md_a = "Ana_Manu_sol_N"*string(N)*"_K"*string(K)*".mat"
        ref_data = matread(sol_md_a)
        x_ref = ref_data["x"]
        rho_ref = ref_data["rho_a"]
        rhou_ref = ref_data["rhou_a"] 
        E_ref = ref_data["E_a"]      
        num_data = matread(sol_md)
        x_num = num_data["x"]
        wJq = num_data["wJq"]
        num_rho_N = num_data["rho"]
        err = 0.0
        for e in axes(x_num, 2)
            for ii in axes(x_num, 1)
                err += wJq[ii, e] * (num_rho_N[ii, e] - rho_ref[ii, e])^2
            end
        end
        @show err
        error_mat[i, j] = sqrt(err)
    end
end

#@show error_mat
e_mat = error_mat[2, :]
println("Error mat is: $e_mat")
r_n1 = (log(error_mat[1, end-1]/error_mat[1, end]))/log(2)
r_n2 = (log(error_mat[2, end-1]/error_mat[2, end]))/log(2)
r_n3 = (log(error_mat[3, end-1]/error_mat[3, end]))/log(2)
r_n5 = (log(error_mat[5, end-1]/error_mat[5, end]))/log(2)
@show error_mat

@show r_n1
@show r_n2
@show r_n3
@show r_n4
@show r_n5

K = [2, 4, 8, 16, 32]

plot(K, error_mat[1, :], label="N=1")
plot!(K, error_mat[2, :], label="N=2")
plot!(K, error_mat[3, :], label="N=3")
plot!(K, error_mat[4, :], label="N=4")
plot!(K, error_mat[5, :], label="N=5")
plot!(xscale=:log10, yscale=:log10, minorgrid=true)
