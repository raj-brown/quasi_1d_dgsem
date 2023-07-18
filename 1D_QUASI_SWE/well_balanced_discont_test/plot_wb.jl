using MAT
using Plots
using LaTeXStrings
plot_font = "Computer Modern"
default(fontfamily=plot_font, framestyle=:box, label=nothing, grid=false)
#scalefontsizes(1.1)

wb_data = matread("wb_test_4_K200.mat");

hu = wb_data["hAu"]./wb_data["A"];
h = wb_data["hA"]./wb_data["A"];
A = wb_data["A"];
b = wb_data["b"];
x = wb_data["x"];
Vp = wb_data["Vp"]
wJq = wb_data["wJq"]


h_act = wb_data["h_0"]
hu_act = wb_data["hu_0"]



## L2 Norm

err = 0.0

for e in axes(x, 2)
    for ii in axes(x, 1)
        err += wJq[ii, e] * (h[ii, e] - h_act[ii, e])^2
    end
end

err = sqrt(err)
@show err

@show "Show L2 err h : %.20f" err


err = 0.0

for e in axes(x, 2)
    for ii in axes(x, 1)
        err += wJq[ii, e] * (hu[ii, e] - hu_act[ii, e])^2
    end
end
err = sqrt(err)

@show "Show L2 err hu : %.20f" err


## L1 Norm
err = 0.0

for e in axes(x, 2)
    for ii in axes(x, 1)
        err += wJq[ii, e] * (h[ii, e] - h_act[ii, e])
    end
end

err = abs(err)
@show err

@show "Show L1 err h : %.20f" err


err = 0.0

for e in axes(x, 2)
    for ii in axes(x, 1)
        err += wJq[ii, e] * (hu[ii, e] - hu_act[ii, e])
    end
end
err = abs(err)


@show "Show L1 err hu : %.20f" err


err = zeros(5, 200)

for e in axes(x, 2)
    for ii in axes(x, 1)
        err[ii, e] = abs(wJq[ii, e] * (h[ii, e] - h_act[ii, e]))
    end
end

err = max(vcat(err)...)
@show err

@show "Show Linfty err h: %.20f" err


err = zeros(5, 200)

for e in axes(x, 2)
    for ii in axes(x, 1)
        err[ii, e]= abs(wJq[ii, e] * (hu[ii, e] - hu_act[ii, e]))
    end
end
err = max(vcat(err)...)
@show "Show Linfty err hu : %.20f" err






xp = Vp * x
b =  Vp * b
h = Vp * h
A = Vp * A


xp = vcat(xp...)
b = vcat(b...)
h = vcat(h...)
A = vcat(A...)


plot(xp, h .+ b, linewidth=3, c=:red, 
legend=:topleft, labels="water surface (h+b)", xtickfontsize=14, ytickfontsize=14)
plot!(xp, b, linewidth=3, c=:blue, 
legend=:topleft, labels="bottom (b)", xtickfontsize=14, ytickfontsize=14)
ylims!((0,1.5))
xlabel!("x")
ylabel!("water surface, bottom")

savefig("top_bottom.png") 


plot(xp, A, linewidth=3, c=:black, xtickfontsize=14, ytickfontsize=14)
#scatter!([0.5], [0.6], mc=:red, ms=5, ma=2.0, legend=:bottomleft, labels=L"$1 - 2 \sigma_0$", xtickfontsize=14, ytickfontsize=14)
yticks!([0.0,0.2,0.4,0.6, 0.8, 1.0])
xlabel!(L"x")
ylabel!(L"a(x)")
ylims!((0,1))

savefig("channel.png") 


# Entropy Data
# ec_data = matread("ec_t.mat")
# N_skip =100

# ### Extract analytical data
# t = vec(ec_data["t"])
# ec = vec(ec_data["ec_val"])
# x = ec_data["x"]
# x = [x[1, :]; x[end, end]]

# t = t[1:N_skip:end ]
# ec = ec[1:N_skip:end ]

# plot(t, ec, linewidth=3, c=:red, legend=:false)
# ylims!((10^-30, 10^0))
# plot!( yscale=:log10, minorgrid=true)

# xlabel!("Time (t)")
# ylabel!("Entropy residuals")
# #savefig("ec.png") 