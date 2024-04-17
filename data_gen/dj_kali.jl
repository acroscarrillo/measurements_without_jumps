include("../src/src.jl")

using DataFrames
using ProgressBars


N = 500
Δ = 1000
χ_x = 0
χ_y = 1
χ_z = 1

Tmax = 2*π
Nt = 100
t_array = Vector(range(0,Tmax,Nt))
dt = Tmax/Nt

σ_x_2,σ_x_3,σ_x_4 =  σ_n_j(4,1,2),σ_n_j(4,1,3),σ_n_j(4,1,4)


# ψ_1_0 = [1, 2*exp(im*π/8)]/√(5)

ψ_1_0 = [1, 0]

n = [ψ_1_0'*σ_1*ψ_1_0 ψ_1_0'*σ_2*ψ_1_0 ψ_1_0'*σ_3*ψ_1_0]

# ψ_1_0 = [1, 0]

# ψ_0 = ψ_init(ψ_1_0)

plus_x = [1, 1]/√(2)
plus_z = [0, 1]
plus_y = [1, im]/√(2)
minus_z = [1, 0]

# ψ_0 = kron( kron( kron(plus_z,plus_x) ,plus_y) ,plus_y)
# ψ_0 = kron( kron( kron(minus_z,plus_x) ,plus_y) ,plus_y)

Bell_PHI_plus = (kron(minus_z,minus_z)+kron(plus_z,plus_z))/sqrt(2)
Bell_PSI_plus = (kron(minus_z,plus_z)+kron(plus_z,minus_z))/sqrt(2)

# ψ_0 = kron( kron(Bell_PHI_plus ,plus_y) ,plus_y)
ψ_0 = kron( kron(Bell_PSI_plus ,plus_y) ,plus_y)

ψ_t = ψ_0
measured_data = zeros(length(t_array)*N, 6)
overlap_array = zeros(length(t_array))
counter = 1
for (j,t) in ProgressBar(enumerate(t_array))
    U_t = exp(-im*t*H(n,Δ,χ_x,χ_y,χ_z))
    for n=1:N
        ψ_t = U_t*ψ_t

        x_2, ψ_t = measure(ψ_t,σ_x_2)
        x_3, ψ_t = measure(ψ_t,σ_x_3)
        x_4, ψ_t = measure(ψ_t,σ_x_4)

        ψ_t = cond_uni(x_2,2)*ψ_t
        ψ_t = cond_uni(x_3,3)*ψ_t
        ψ_t = cond_uni(x_4,4)*ψ_t

        overlap = norm(ψ_t'*ψ_0)

        measured_data[counter,:] .= overlap,x_2,x_3,x_4,n,t

        # display((ψ_t'*σ_n_j(4,2,2)*ψ_t,ψ_t'*σ_n_j(4,2,3)*ψ_t,ψ_t'*σ_n_j(4,2,4)*ψ_t))
        counter += 1
    end
    
end

data_df = DataFrame(measured_data, ["overlap","x_2","x_3","x_4","n","t"])

avg_x_2_t = zeros(length(t_array))
avg_x_3_t = zeros(length(t_array))
avg_x_4_t = zeros(length(t_array))
overlap_t = zeros(length(t_array))
for j=1:length(t_array)
    df_temp = filter(row ->  row.t == t_array[j], data_df)
    avg_x_2_t[j] = mean(df_temp.x_2)
    avg_x_3_t[j] = mean(df_temp.x_3)
    avg_x_4_t[j] = mean(df_temp.x_4)
    overlap_t[j] = mean(df_temp.overlap)
end

plot(t_array, avg_x_2_t,ylim = (-1.1,1.1))
scatter!(t_array, avg_x_2_t,ylim = (-1.1,1.1))
plot!(t_array, avg_x_3_t)
scatter!(t_array, avg_x_3_t)
plot!(t_array, avg_x_4_t)
scatter!(t_array, avg_x_4_t)
plot!(t_array, overlap_t)
scatter!(t_array, overlap_t)




