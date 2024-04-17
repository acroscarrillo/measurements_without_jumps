include("../src/src.jl")

using DataFrames
using ProgressBars
using Plots
using LaTeXStrings

σ_x_2,σ_x_3,σ_x_4 =  σ_n_j(4,1,2),σ_n_j(4,1,3),σ_n_j(4,1,4)

#########################
# Simulation parameters #
#########################
N_tot_runs = 10000

N = 100
Δ_array = Vector(range(20,200,length=10))
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
final_fidelity_array = zeros(N_tot_runs,length(Δ_array))
pbar = ProgressBar(total=length(final_fidelity_array))
for (k,Δ) in enumerate(Δ_array)
    for m=1:N_tot_runs
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
        final_fidelity_array[m,k] = norm(ψ_t'*ψ_0)^2   
        update(pbar)
    end
end

P_threshold = 0.99
success_array =  zeros(length(Δ_array))
for (k,Δ) in enumerate(Δ_array)
    temp = final_fidelity_array[:,k] 
    success_array[k] = length(temp[ temp .> P_threshold ])/ N_tot_runs
end


width_px= 340.39020340390204  # or  width_inches*DPI
scatter(Δ_array,success_array,size = (width_px,width_px*0.6),legend=false,ylabel=L"\textrm{Success \ \ probability}",xlab=L"\Delta/\chi",xtickfontsize=8,ytickfontsize=8,guidefont=font(10),dpi=600,widen=false,tickdirection=:out,right_margin = 2Plots.mm,bottom_margin = 0Plots.mm,fontfamilys = "Times New Roman",legendfontsize=9,ms=3,xlim=(-10 + Δ_array[1],10 + Δ_array[end]),ylim=(0.6,1.05),title=L"N = "*"$N,   "*L"N_t="*"$Nt",alpha=0.75,titlefontsize=10)

# savefig("figs/success_prob_vs_Delta.png")
# savefig("figs/success_prob_vs_Delta.pdf")