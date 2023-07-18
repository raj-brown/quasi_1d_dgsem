using MAT
using Plots
using LaTeXStrings
using CSV, DataFrames
plot_font = "Computer Modern"
default(fontfamily=plot_font, framestyle=:box, label=nothing, grid=false)
#scalefontsizes(1.1)

wb_data = matread("Sol_N2_K200.mat");

hu = wb_data["hAu"]./wb_data["A"];
h = wb_data["hA"]./wb_data["A"];
A = wb_data["A"];
b = wb_data["b"];
x = wb_data["x"];
Vp = wb_data["Vp"]
wJq = wb_data["wJq"]
u = hu./h
H = A.*h
g = 9.81
c = sqrt.(g * H./ A)

f = abs.(u)./c


x = Vp * x
b =  Vp * b
h = Vp * h
A = Vp * A

f = Vp * f


xp = vcat(x...)
b = vcat(b...)
h = vcat(h...)
A = vcat(A...)
f = vcat(f...)


h_vq_d = CSV.read("fig_23.csv", DataFrame, header=false, missingstring="NA", delim=',')
f_vq_d = CSV.read("fig_24.csv", DataFrame, header=false, missingstring="NA", delim=',')


x_vq = h_vq_d[!,"Column1"]
h_vq = h_vq_d[!,"Column2"]


N_skip = 8
plot(xp[1:N_skip:end], h[1:N_skip:end], linewidth=3.0, c=:black, legend=:topleft, labels="ESDG Method", xtickfontsize=14, ytickfontsize=14)
scatter!(x_vq[1:2:end], h_vq[1:2:end], markersize=7, markercolor=:blue, shape=:circle, markeralpha=0.3, legend=:bottomleft, labels="Vázquez-Céndon (1999)", xtickfontsize=14, ytickfontsize=14)
yticks!([0.0,0.5,1.0,1.5, 2.0, 2.5, 3.0])
ylims!((0, 2.1))
xlabel!("x")
ylabel!(L"h")
savefig("h.png") 


x_vq = f_vq_d[!,"Column1"]
f_vq = f_vq_d[!,"Column2"]


plot(xp[1:N_skip:end], f[1:N_skip:end], linewidth=3.0, c=:black, legend=:topleft, labels="ESDG Method", xtickfontsize=14, ytickfontsize=14)
scatter!(x_vq[1:2:end], f_vq[1:2:end], markersize=7, markercolor=:blue, shape=:circle, markeralpha=0.3, legend=:topleft, labels="Vázquez-Céndon (1999)", xtickfontsize=14, ytickfontsize=14)
yticks!([0.0,0.5,1.0,1.5, 2.0, 2.5, 3.0])
ylims!((0, 2.1))
xlabel!("x")
ylabel!(L"Froude number: $\frac{|u|}{\sqrt{gh}}$")
savefig("f.png") 


