using LinearAlgebra
using StatsBase
# struct Circuit
#     ...

global σ_0 = [1 0; 0 1]
global σ_1 = [0 1; 1 0]
global σ_2 = [0 -im; im 0]
global σ_3 = [1 0; 0 -1]

function σ_n(n)
    if n == 0
        return σ_0
    elseif n == 1
        return σ_1
    elseif n == 2
        return σ_2
    elseif n == 3
        return σ_3
    else
        return 0 
    end
end

function σ_n_j(N,n,j)
    op_list = []
    for i in range(1,N)
        if i == j
            push!(op_list,σ_n(n))
        else
            push!(op_list,σ_n(0))
        end
    end
    return tensor_list(op_list)
end

function tensor_list(list_of_ops)
    prod = 1
    for op in list_of_ops
        prod = kron(prod,op)
    end
    return prod
end

function measure(ψ,Op,deg_thresh=1e-3)
    lambs, vec = eigen(Op)
    prob_list = (norm.(vec'*ψ)).^2
    outcomes = []
    for i in range(1,length(lambs))
        push!(outcomes, (lambs[i],vec[:,i]*vec[:,i]'))
    end
    λ_obs, P_obs =  sample(outcomes, Weights(prob_list))

    P = zeros(length(ψ),length(ψ))
    for j=1:length(ψ)
        λ = lambs[j]
        δλ = abs(λ - λ_obs)
        if δλ < deg_thresh
            P += vec[:,j]*vec[:,j]'
        end
    end

    ψ_collapsed = P*ψ
    return λ_obs, ψ_collapsed/norm(ψ_collapsed)
end


function H(n,Δ,χ_x,χ_y,χ_z)
    n = n/norm(n)
    H_0 = Δ * ( n[1]*σ_n_j(4,1,1) + n[2]*σ_n_j(4,2,1) + n[3]*σ_n_j(4,3,1) )
    H_xz = χ_x * σ_n_j(4,1,1)*σ_n_j(4,3,2) 
    H_yz = χ_y * σ_n_j(4,2,1)*σ_n_j(4,3,3) 
    H_zz = χ_z * σ_n_j(4,3,1)*σ_n_j(4,3,4) 
    return H_0 + H_xz + H_yz + H_zz
end
