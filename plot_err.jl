using LinearAlgebra
using Plots
using MAT
using Revise
using Interpolations

error_mat = zeros(5, 5)

for i=1:5
    for j =1:5
        N = i
        K = 2^j
        ref_data = matread("/Users/raj/WORK_RAJ/DG_SEM/code/Ref_Sol_N_5_K_2000.mat")
        @show keys(ref_data)
        x_ref = ref_data["x"]
        rho_ref = ref_data["rho"]
        rhou_ref = ref_data["rhou"] 
        E_ref = ref_data["E"] 
        plot(x_ref[1, :], rho_ref[1, :], leg=false)
        gui()
       
        x_ref = [round(x_ref[1, 1]); x_ref[end, :]]
        rho_ref = [rho_ref[1, 1]; rho_ref[end, :]]
        fn_rho_ref = linear_interpolation(x_ref, rho_ref)
        y = fn_rho_ref.(x_ref)
        plot!(x_ref, y, ls=:dot)
        gui()
        f_name = "sol_N"*string(N)*"_K"*string(K)*".mat"
        num_data = matread(f_name)
        x_num = num_data["x"]

        wJq = num_data["wJq"]
        num_rho_N = num_data["rho"]
        err = 0.0
        for e in axes(x_num, 2)
            for ii in axes(x_num, 1)
                err += wJq[ii, e] * (num_rho_N[ii, e] - fn_rho_ref(x_num[ii, e]))^2
            end
        end
        error_mat[i, j] = sqrt(err)

        # x_num = [round(x_num[1, 1]); x_num[end, :]]

        # ref_rho_N = fn_rho_ref(x_num)
        # num_rho_N = num_data["rho"]
        # num_rho_N = [num_rho_N[1, 1]; num_rho_N[end, :]]


        # @infiltrate
        # err_2 = sqrt(sum(wJq .* (ref_rho_N - num_rho_N).^2)) #norm(ref_rho_N - num_rho_N) / norm(ref_rho_N)
        # @show N, err_2
        # error_mat[i, j] = err_2
    end
end

#@show error_mat
e_mat = error_mat[2, :]
println("Error mat is: $e_mat")
r_n1 = (log(error_mat[1, end-1]/error_mat[1, end]))/log(2)
r_n2 = (log(error_mat[2, end-1]/error_mat[2, end]))/log(2)
r_n3 = (log(error_mat[3, end-1]/error_mat[3, end]))/log(2)
r_n4 = (log(error_mat[4, end-2]/error_mat[4, end-1]))/log(2)
r_n5 = (log(error_mat[5, end-2]/error_mat[5, end-1]))/log(2)
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
gui()