include("../src/src.jl")

using DataFrames
using ProgressBars
using Plots

σ_x_2,σ_x_3,σ_x_4 =  σ_n_j(4,1,2),σ_n_j(4,1,3),σ_n_j(4,1,4)

#########################
# Simulation parameters #
#########################

N = 5000
Δ = 200
χ_x = 1
χ_y = 1
χ_z = 1

Tmax = 2*π
Nt = 100
t_array = Vector(range(0,Tmax,Nt))
dt = Tmax/Nt


#################
# Initial state #
#################

ψ_1_0 = [1, 2*exp(im*π/8)]/√(5)

# ψ_1_0 = [1, 0]

n = [ψ_1_0'*σ_1*ψ_1_0 ψ_1_0'*σ_2*ψ_1_0 ψ_1_0'*σ_3*ψ_1_0]

# ψ_1_0 = [1, 0]

ψ_0 = ψ_init(ψ_1_0) # takes vector and tensors it with y_plus = [1, im]/√(2)


#################
# Initial state #
#################
ψ_t = ψ_0
measured_data = zeros(length(t_array)*N, 6)
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

overlap_t_avg = zeros(length(t_array))
overlap_t_err = zeros(length(t_array))
for j=1:length(t_array)
    df_temp = filter(row ->  row.t == t_array[j], data_df)
    overlap_t_avg[j] = mean(df_temp.overlap)
end

plot(t_array, overlap_t_avg, yscale=:log10)
scatter!(t_array, overlap_t_avg)
