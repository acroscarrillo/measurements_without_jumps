include("../src/src.jl")

using DataFrames
using ProgressBars
using Plots
using LaTeXStrings

σ_x_2,σ_x_3,σ_x_4 =  σ_n_j(4,1,2),σ_n_j(4,1,3),σ_n_j(4,1,4)

#########################
# Simulation parameters #
#########################
N_tot_runs = 1000

N = 100
Δ = 100
Nt = 4

χ_x = 1
χ_y = 1
χ_z = 1

Tmax = π
t_array = Vector(range(0,Tmax,Nt))
dt = Tmax/Nt

#################
# Initial state #
#################

ψ_1_0 = [1, 2*exp(im*π/8)]/√(5)

# ψ_1_0 = [1, 0]

n = real.([ψ_1_0'*σ_1*ψ_1_0 ψ_1_0'*σ_2*ψ_1_0 ψ_1_0'*σ_3*ψ_1_0])

# ψ_1_0 = [1, 0]

ψ_0 = ψ_init(ψ_1_0) # takes vector and tensors it with y_plus = [1, im]/√(2)

#################
# Initial state #
#################
final_fidelity = zeros(N_tot_runs)
for m=ProgressBar(1:N_tot_runs)
    ψ_t = ψ_0
    for (j,t) in enumerate(t_array)
        U_t = exp(-im*(0-t)*H(n,Δ,χ_x,χ_y,χ_z))
        # measurements:
        for n=1:N
            ψ_t = U_t*ψ_t

            x_2, ψ_t = measure(ψ_t,σ_x_2)
            x_3, ψ_t = measure(ψ_t,σ_x_3)
            x_4, ψ_t = measure(ψ_t,σ_x_4)

            ψ_t = cond_uni(x_2,2)*ψ_t
            ψ_t = cond_uni(x_3,3)*ψ_t
            ψ_t = cond_uni(x_4,4)*ψ_t

        end
    end
    final_fidelity[m] = norm(ψ_t'*ψ_0)^2    
end

P_threshold = 0.5
prob = round(length(final_fidelity[ final_fidelity .> P_threshold ])/ N_tot_runs, sigdigits=3)

width_px= 340.39020340390204  # or  width_inches*DPI
scatter(final_fidelity,size = (width_px,width_px*0.6),label=false,ylabel=L"\mathcal{F}_f",xlab=L"\textrm{Trails}",xtickfontsize=8,ytickfontsize=8,guidefont=font(10),dpi=600,widen=false,tickdirection=:out,right_margin = 2Plots.mm,bottom_margin = 0Plots.mm,fontfamilys = "Times New Roman",legendfontsize=9,ms=3,ylim=(-0.1,1.1),title=L"N = "*"$N,   "*L"\Delta/\chi="*"$Δ,   "*L"\mathbb{P}[\mathcal{F}_f >"*"$P_threshold "*L"] = "*"$prob",alpha=0.75,titlefontsize=10)

hline!([P_threshold],lw=3, c=:red,label=L"\mathcal{F}_f = "*"$P_threshold",ls=:dash,legend=:left)

# savefig("figs/fail_prob_N_100_Delta_50.png")
# savefig("figs/fail_prob_N_100_Delta_50.pdf")