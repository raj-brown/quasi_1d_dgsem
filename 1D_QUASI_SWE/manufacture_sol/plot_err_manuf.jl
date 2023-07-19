using LinearAlgebra
using Plots
using MAT

error_mat_h = zeros(5, 5)
error_mat_hu = zeros(5, 5)




N_arr = [1, 2, 3, 4, 5]
K_arr = [16, 32, 64, 128, 256]

for i=eachindex(N_arr)
    for j =eachindex(K_arr)
        
        N = N_arr[i]
        K = K_arr[j]
        @show N, K
        sol_md = "Num_Manu_sol_N"*string(N)*"_K"*string(K)*".mat"
        sol_md_a = "Ana_Manu_sol_N"*string(N)*"_K"*string(K)*".mat"
        
        #####
        ref_data = matread(sol_md_a)
        x_ref = ref_data["x"]
        ha_ref = ref_data["hA_a"]
        hua_ref = ref_data["hAu_a"]
        A_a = ref_data["A_a"] 
        h_ref = ha_ref./A_a
        hu_ref = hua_ref./A_a

      
        ######
        num_data = matread(sol_md)
        @show keys(num_data)
        x_num = num_data["x"]
        wJq = num_data["wJq"]
        num_ha_N = num_data["hA"]
        num_hua_N = num_data["hAu"]
        A_num = num_data["A"] 

        h_num =  num_ha_N./A_num
        hu_num = num_hua_N./A_num

        err_1 = 0.0
        err_2 = 0.0
        for e in axes(x_num, 2)
            for ii in axes(x_num, 1)
                err_1 += wJq[ii, e] * (h_num[ii, e] - h_ref[ii, e])^2
                err_2 += wJq[ii, e] * (hu_num[ii, e] - hu_ref[ii, e])^2
            end
        end
        error_mat_h[i, j] = sqrt(err_1)
        error_mat_hu[i, j] = sqrt(err_2)
    end
end

error_mat = sqrt.(error_mat_h.^2 + error_mat_hu.^2)

r_n1 = (log(error_mat[1, end-1]/error_mat[1, end]))/log(2)
r_n2 = (log(error_mat[2, end-1]/error_mat[2, end]))/log(2)
r_n3 = (log(error_mat[3, end-1]/error_mat[3, end]))/log(2)
r_n4 = (log(error_mat[4, end-1]/error_mat[4, end]))/log(2)
r_n5 = (log(error_mat[5, end-1]/error_mat[5, end]))/log(2)

@show r_n1
@show r_n2
@show r_n3
@show r_n4
@show r_n5

plot(K_arr, error_mat[1, :], label="N=1");
plot!(K_arr, error_mat[2, :], label="N=2");
plot!(K_arr, error_mat[3, :], label="N=3");
plot!(K_arr, error_mat[4, :], label="N=4");
plot!(K_arr, error_mat[5, :], label="N=5");
plot!(xscale=:log10, yscale=:log10, minorgrid=true);

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


