using Plots
using MAT
using Interpolations

error_mat_rho = zeros(5, 5)
error_mat_rhou = zeros(5, 5)
error_mat_E = zeros(5, 5)



for i=1:5
    for j =1:5
        N = i
        K = 2^j
        ref_data = matread("Rev_Ref_Sol_N_3_K_8000.mat")
        @show keys(ref_data)
        x_ref = ref_data["x"]
        rho_ref = ref_data["rho"]
        rhou_ref = ref_data["rhou"] 
        E_ref = ref_data["E"] 
        plot(x_ref[1, :], rho_ref[1, :], leg=false)
        gui()
       
        x_ref = [round(x_ref[1, 1]); x_ref[end, :]]
        rho_ref = [rho_ref[1, 1]; rho_ref[end, :]]
        rhou_ref = [rhou_ref[1, 1]; rhou_ref[end, :]]
        E_ref = [E_ref[1, 1]; E_ref[end, :]]

        fn_rho_ref = linear_interpolation(x_ref, rho_ref)
        fn_rhou_ref = linear_interpolation(x_ref, rhou_ref)
        fn_E_ref = linear_interpolation(x_ref, E_ref)

        f_name = "sol_N"*string(N)*"_K"*string(K)*".mat"
        num_data = matread(f_name)
        x_num = num_data["x"]

        wJq = num_data["wJq"]
        num_rho_N = num_data["rho"]
        num_rhou_N = num_data["rhou"]
        num_E_N = num_data["E"]

        err_1 = 0.0
        err_2 = 0.0
        err_3 = 0.0
        for e in axes(x_num, 2)
            for ii in axes(x_num, 1)
                err_1 += wJq[ii, e] * (num_rho_N[ii, e] - fn_rho_ref(x_num[ii, e]))^2
                err_2 += wJq[ii, e] * (num_rhou_N[ii, e] - fn_rhou_ref(x_num[ii, e]))^2
                err_3 += wJq[ii, e] * (num_E_N[ii, e] - fn_E_ref(x_num[ii, e]))^2

            end
        end
        error_mat_rho[i, j] = sqrt(err_1)
        error_mat_rhou[i, j] = sqrt(err_2)
        error_mat_E[i, j] = sqrt(err_3)


    end
end
error_mat = sqrt.(error_mat_rho.^2 + error_mat_rhou.^2 + error_mat_E.^2)

#@show error_mat
println("Error mat is: $error_mat")
r_n1 = (log(error_mat[1, end-1]/error_mat[1, end]))/log(2)
r_n2 = (log(error_mat[2, end-1]/error_mat[2, end]))/log(2)
r_n3 = (log(error_mat[3, end-1]/error_mat[3, end]))/log(2)
r_n4 = (log(error_mat[4, end-1]/error_mat[4, end-1]))/log(2)
r_n5 = (log(error_mat[5, end-1]/error_mat[5, end-1]))/log(2)

@show r_n1
@show r_n2
@show r_n3
@show r_n4
@show r_n5

K = [2, 4, 8, 16, 32]

r_N1 = log.(error_mat[1, 1:end-1]./error_mat[1, 2:end])./log(2) ;
r_N2 = log.(error_mat[2, 1:end-1]./error_mat[2, 2:end])./log(2) ;
r_N3 = log.(error_mat[3, 1:end-1]./error_mat[3, 2:end])./log(2) ;
r_N4 = log.(error_mat[4, 1:end-1]./error_mat[4, 2:end])./log(2) ;
r_N5 = log.(error_mat[5, 1:end-1]./error_mat[5, 2:end])./log(2) ;


@show r_N1
@show r_N2
@show r_N3
@show r_N4
@show r_N5



plot(K, error_mat[1, :], label="N=1")
plot!(K, error_mat[2, :], label="N=2")
plot!(K, error_mat[3, :], label="N=3")
plot!(K, error_mat[4, :], label="N=4")
plot!(K, error_mat[5, :], label="N=5")
plot!(xscale=:log10, yscale=:log10, minorgrid=true)

