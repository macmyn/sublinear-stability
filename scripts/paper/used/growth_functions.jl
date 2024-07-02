using PyPlot, Revise, DrWatson, OMEinsum, Random,LaTeXStrings, ForwardDiff

include(srcdir("NonlinearStability.jl"))

## USED IN PAPER showing behaviour of growth functions near 0

plt.rc("text", usetex=true)
rc("font", family="serif", weight="normal", size="18")
rc("axes", labelsize="10")
rc("xtick", labelsize="10")
rc("ytick", labelsize="10")

xs = 0:0.01:2

ys_lv_pos(x) =  x*(1-x^(0.5))
ys_lv_neg(x) =  x*(-1+x^(-0.5))

d_y_pos(x) = ForwardDiff.derivative(ys_lv_pos,x)
d_y_neg(x) = ForwardDiff.derivative(ys_lv_neg,x)

fig, ax = plt.subplots(figsize = (2,1.6))

axhline(0,color="k",linewidth=0.5) # x = 0
axvline(0,color="k",linewidth=0.5) # y = 0

cm = get_cmap(:tab10)
colorrange = (0:9)./10

ax.plot(xs,ys_lv_pos.(xs),label=L"\alpha=0.5",color=cm(colorrange[1]))
ax.plot(xs,ys_lv_neg.(xs),label=L"\alpha=-0.5",color=cm(colorrange[2]))
ax.plot(xs,d_y_pos.(xs),color=cm(colorrange[1]),linestyle="--",alpha=0.5)
ax.plot(xs,d_y_neg.(xs),color=cm(colorrange[2]),linestyle="--",alpha=0.5)

ax.set_ylim(-0.5,1.5)
ax.set_xlim(0,)

ax.set_xlabel(L"N")
ax.set_ylabel(L"\frac{\mathrm{d}N}{\mathrm{dt}}")

ax.legend(fontsize=10)

savefig(plotsdir("alpha_pos_neg_deriv.pdf"),bbox_inches="tight")