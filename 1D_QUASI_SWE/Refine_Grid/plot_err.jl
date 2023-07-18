using LinearAlgebra
using Plots
using MAT
using Interpolations

error_mat_h = zeros(5, 6)
error_mat_hu = zeros(5, 6)


K_arr = [2, 4, 8, 16, 32, 64]
#K_arr = [10, 20, 40, 80, 160]

N_arr = [1, 2, 3, 4, 5]

for i=1:5
    for j =1:6
        N = N_arr[i]
        K = K_arr[j]
        ref_data = matread("/Users/raj/WORK_RAJ/DG_SEM/code/1D_quasi_swe/Refine_grid/Ref_Sol_N3_K8000.mat")
        @show keys(ref_data)
        x_ref = ref_data["x"]
        hA_ref = ref_data["hA"]
        hAu_ref = ref_data["hAu"] 
        A = ref_data["A"]
        h_ref = hA_ref./A 
        hu_ref = hAu_ref./A
        plot(x_ref[1, :], h_ref[1, :], leg=false)
        #gui()
       
        x_ref = [round(x_ref[1, 1]); x_ref[end, :]]
        #@show "x_ref is", x_ref[1]
        h_ref = [h_ref[1, 1]; h_ref[end, :]]
        hu_ref = [hu_ref[1, 1]; hu_ref[end, :]]

        fn_h_ref = linear_interpolation(x_ref, h_ref)
        fn_hu_ref = linear_interpolation(x_ref, hu_ref)
        #fn_h_ref = CubicSplineInterpolation(x_ref, h_ref)
        #fn_hu_ref = CubicSplineInterpolation(x_ref, hu_ref)

        f_name = "Ref_Sol_N"*string(N)*"_K"*string(K)*".mat"
        num_data = matread(f_name)
        x_num = num_data["x"]
        #x_num = [round(x_num[1, 1]); x_num[end, :]]
        x_num[1,1] = 0.0
        wJq = num_data["wJq"]
        num_h_N  = num_data["hA"]./num_data["A"]
        num_hu_N = num_data["hAu"]./num_data["A"]

        err_1 = 0.0
        err_2 = 0.0
        for e in axes(x_num, 2)
            for ii in axes(x_num, 1)
                err_1 += wJq[ii, e] * (num_h_N[ii, e] - fn_h_ref(x_num[ii, e]))^2
                err_2 += wJq[ii, e] * (num_hu_N[ii, e] - fn_hu_ref(x_num[ii, e]))^2
            end
        end
        error_mat_h[i, j] = sqrt(err_1)
        error_mat_hu[i, j] = sqrt(err_2)
    end
end
error_mat = sqrt.(error_mat_h.^2 + error_mat_hu.^2)

#@show error_mat
#println("Error mat is: $e_mat")
r_n1 = (log(error_mat[1, end-1]/error_mat[1, end]))/log(2)
r_n2 = (log(error_mat[2, end-1]/error_mat[2, end]))/log(2)
r_n3 = (log(error_mat[3, end-1]/error_mat[3, end]))/log(2)
r_n4 = (log(error_mat[4, end-1]/error_mat[4, end]))/log(2)
r_n5 = (log(error_mat[5, end-1]/error_mat[5, end]))/log(2)
@show error_mat

@show r_n1
@show r_n2
@show r_n3
@show r_n4
@show r_n5

K = [2, 4, 8, 16, 32, 64]

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

