"""
ψ_init(ψ_1)

Helper function to append to the given initial state in the first qubit 3 qubits in the |+y> state. In other words, returns:

    |ψ_0> = |ψ_1> ⨂ |+y> ⨂ |+y> ⨂ |+y> 

where

    |+y> = [1, im]/√(2)

# Examples
```julia-repl
julia> round.(  ψ_init([1,0])', sigdigits=3)
1×16 Matrix{ComplexF64}:
    0.354-0.0im  0.0-0.354im  0.0-0.354im  -0.354-0.0im  0.0-0.354im  -0.354-0.0im  -0.354-0.0im  -0.0+0.354im  0.0-0.0im  0.0-0.0im  0.0-0.0im  0.0-0.0im  0.0-0.0im  0.0-0.0im  0.0-0.0im  0.0-0.0im
````
"""
function ψ_init(ψ_1)
    y_plus = [1, im]/√(2)
    temp = ψ_1
    for _=1:3
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