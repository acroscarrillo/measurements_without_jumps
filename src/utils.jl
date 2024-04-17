
function ψ_init(ψ_1)
    y_plus = [1, im]/√(2)
    temp = ψ_1
    for i=1:3
        temp = kron(temp,y_plus)
    end
    return temp
end

function cond_uni(measured_outcome,j)
    outcome = sign(measured_outcome)
    op_list = []
    if outcome == 1
        U = [1 0; 0 im]
    elseif outcome == -1
        U = [1 0; 0 -im]
    else 
        throw(ErrorException("measured $outcome, which is not 1, -1"))
    end

    for n=1:4
        if n == j
            push!(op_list, U)
        else
            push!(op_list, I(2))
        end
    end
    return tensor_list(op_list)
end