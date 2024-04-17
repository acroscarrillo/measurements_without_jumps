include("../src/src.jl")

using DataFrames
using ProgressBars
using Plots
using LaTeXStrings

σ_x_2,σ_x_3,σ_x_4 =  σ_n_j(4,1,2),σ_n_j(4,1,3),σ_n_j(4,1,4)

#########################
# Simulation parameters #
#########################

N = 100
Δ = 50
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
ψ_t = ψ_0
measured_data = zeros(length(t_array)*N, 6)
Fidelity_array = zeros(length(t_array))
counter = 1
for (j,t) in ProgressBar(enumerate(t_array))
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

        Fidelity = norm(ψ_t'*ψ_0)^2

        measured_data[counter,:] .= Fidelity,x_2,x_3,x_4,n,t

        # display((ψ_t'*σ_n_j(4,2,2)*ψ_t,ψ_t'*σ_n_j(4,2,3)*ψ_t,ψ_t'*σ_n_j(4,2,4)*ψ_t))
        counter += 1
    end
    
end

data_df = DataFrame(measured_data, ["overlap","x_2","x_3","x_4","n","t"])

avg_x_2_t_avg = zeros(length(t_array))
avg_x_3_t_avg = zeros(length(t_array))
avg_x_4_t_avg = zeros(length(t_array))
avg_x_2_t_err = zeros(length(t_array))
avg_x_3_t_err = zeros(length(t_array))
avg_x_4_t_err = zeros(length(t_array))

overlap_t = zeros(length(t_array))
for j=1:length(t_array)
    df_temp = filter(row ->  row.t == t_array[j], data_df)
    avg_x_2_t_avg[j] = mean(df_temp.x_2)
    avg_x_3_t_avg[j] = mean(df_temp.x_3)
    avg_x_4_t_avg[j] = mean(df_temp.x_4)
    avg_x_2_t_err[j] = std(df_temp.x_2)/√(N)
    avg_x_3_t_err[j] = std(df_temp.x_3)/√(N)
    avg_x_4_t_err[j] = std(df_temp.x_4)/√(N)

    overlap_t[j] = mean(df_temp.overlap)
end



t_array_normalised = t_array/(π*1) #careful, xi is hard coded!!!

# fidelity 
width_px= 340.39020340390204  # or  width_inches*DPI

# theory lines
t_theory_array = Vector(range(0,1,1000))

plot(t_theory_array, sin.(2*π*t_theory_array.*n[1]),c = :orange,lw=3,label=false)
plot!(t_theory_array, sin.(2*π*t_theory_array.*n[2]),c = :navy,lw=3,label=false)
plot!(t_theory_array, sin.(2*π*t_theory_array.*n[3]),c = :magenta,lw=3,label=false)


scatter!(t_array_normalised, avg_x_2_t_avg,yerr=avg_x_2_t_err,ylim = (-1.1,1.1),c = :orange,label=L"\langle \sigma^{(1)}_x \rangle",ms=3)
scatter!(t_array_normalised, avg_x_3_t_avg,yerr=avg_x_3_t_err,c = :navy,label=L"\langle \sigma^{(2)}_x \rangle",ms=3)
scatter!(t_array_normalised, avg_x_4_t_avg,yerr=avg_x_4_t_err,c = :magenta,label=L"\langle \sigma^{(3)}_x \rangle",ms=3)


scatter!(t_array_normalised, overlap_t,size = (width_px,width_px*0.6),ylim=(-1.1,1.1),xtickfontsize=8,ytickfontsize=8,guidefont=font(10),dpi=600,widen=false,tickdirection=:out,right_margin = 0Plots.mm,left_margin = 0Plots.mm,bottom_margin = 0Plots.mm,fontfamilys = "Times New Roman",xlim=(-0.05 + t_array_normalised[1],t_array_normalised[end] + 0.05), legendfontsize=9,ylabel=L"\textrm{Observables}",xlabel=L"t \quad (\pi /\chi)",label= L"\overline{\mathcal{F}}^{(0)}",title=L"N = "*"$N,   "*L"N_t="*"$Nt,    "*L"\Delta/\chi="*"$Δ",titlefontsize=10,background_color_legend=:transparent,ms=3,legend = false,m=:xcross,color=:black,markerstrokewidth=3,legend_columns=-1)

savefig("figs/observables_vs_t_N_100_Delta_50.png")
savefig("figs/observables_vs_t_N_100_Delta_50.pdf")