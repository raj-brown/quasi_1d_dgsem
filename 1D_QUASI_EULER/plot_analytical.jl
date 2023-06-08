using MAT
#using MakiePublication
#using CairoMakie

using Plots
using LaTeXStrings

plot_font = "Computer Modern"
default(fontfamily=plot_font, framestyle=:box, label=nothing, grid=false)
scalefontsizes(1.3)



# function myplot(x_a, x_n, M_a, M_n, lab_flow; figure_padding=(2,6,1,6))
#     fig = Figure(figure_padding=figure_padding)
#     ax = Axis(fig, xlabel=L"x~(\text{m})", ylabel=lab_flow, xlabelsize=30, ylabelsize=30,
#     xminorticksvisible=true)


#     # N_skip = 5;
#     # x_n = x_n[1:N_skip:end]
#     # M_n = M_n[1:N_skip:end]


#     s1 = scatter!(ax, x_a, M_a, color="red", markersize=14)
#     l1 = lines!(ax, x_n, M_n)

#     axislegend(ax, [[s1], [l1]], [" ESDG ", "Analytical "],  position=:rt, padding=(0,0,0,0))
    


    
#     fig[1,1] = ax
    
#     return fig

# end

# Subsonic case
ana_data = matread("/Users/raj/WORK_RAJ/DG_SEM/DAVID_MATERIALS/subsonic_analytic.mat");
num_data = matread("/Users/raj/WORK_RAJ/DG_SEM/code/DATA_SUBSONIC/sol_subsonic.mat");

# Trans-sonic case
# ana_data = matread("/Users/raj/WORK_RAJ/DG_SEM/DAVID_MATERIALS/transonic_analytic.mat");
# num_data = matread("/Users/raj/WORK_RAJ/DG_SEM/code/DATA_SUBSONIC/trans_subsonic.mat");



### Extract analytical data
x_a = ana_data["x"];
M_a = ana_data["M"];
p_a = ana_data["p"];

### Vectorized analytical Data
x_a = vec(x_a)
M_a = vec(M_a)
p_a = vec(p_a)

### Extract numerical data
x_n = num_data["x"]
M_n = num_data["Ma"]
p_n = num_data["pressure"]

### Vectorized numerical Data
# x_n = vec(x_n)
# M_n = vec(M_n)
# p_n = vec(p_n)

# plot_sub_sonic_ma() = myplot(x_a, x_n, M_a, M_n, L"Ma" ; figure_padding=(2, 16, 2, 8))
# ma = plot_sub_sonic_ma()
# save("Mach_Trans_Sonic.png", ma, px_per_unit=4)


# plot_sub_sonic_ρ() = myplot(x_a, x_n, p_a, p_n, L"p" ; figure_padding=(2, 16, 2, 8))
# ρ = plot_sub_sonic_ρ()
# save("Pressure_Trans_Sonic.png", ρ, px_per_unit=4)

N_skip = 50

x_a = x_a[1:N_skip:end]
M_a = M_a[1:N_skip:end]
p_a = p_a[1:N_skip:end]


plot(x_n, M_n, linewidth=3, c=:black, legend=:false)
plot(x_n[1:end], M_n[1:end], linewidth=3, c=:black, 
legend=:topleft, labels="ESDG", xtickfontsize=14, ytickfontsize=14)
xlabel!("x (m)")
ylabel!("Ma")


scatter!(x_a, M_a, label="Analytical", mc=:red, ms=5, ma=2.0)
savefig("ss_ma.png") 



plot(x_n, p_n, linewidth=3, c=:black, legend=:false)
plot(x_n[1:end], p_n[1:end], linewidth=3, c=:black, 
legend=:bottomleft, labels="ESDG", xtickfontsize=14, ytickfontsize=14)
xlabel!("x (m)")
ylabel!("p")
#ylim!()


scatter!(x_a, p_a, label="Analytical", mc=:red, ms=5, ma=2.0)
savefig("ss_rho.png") 
