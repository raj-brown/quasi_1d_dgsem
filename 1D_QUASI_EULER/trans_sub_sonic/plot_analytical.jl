using MAT
#using MakiePublication
#using CairoMakie

using Plots
using LaTeXStrings

plot_font = "Computer Modern"
default(fontfamily=plot_font, framestyle=:box, label=nothing, grid=false)
scalefontsizes(1.3)



# Subsonic case
ana_data = matread("subsonic_analytic.mat");
num_data = matread("NUM_DATA/sol_subsonic.mat");

# Trans-sonic case: 
# To generate the result for transonic please comment out the code block: Line: 20-21
# ana_data = matread("transonic_analytic.mat");
# num_data = matread("NUM_DATA/trans_subsonic.mat");



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
savefig("ts_ma.png") 

plot(x_n, p_n, linewidth=3, c=:black, legend=:false)
plot(x_n[1:end], p_n[1:end], linewidth=3, c=:black, 
legend=:bottomleft, labels="ESDG", xtickfontsize=14, ytickfontsize=14)
xlabel!("x (m)")
ylabel!("p")
#ylim!()


scatter!(x_a, p_a, label="Analytical", mc=:red, ms=5, ma=2.0)
savefig("ts_p.png") 
