abstract type AbstractQuantumState{S<:ElementaryQuantumSystem} end

struct GenericQuantumState{K,S,T<:TensorMap{S,K,0},L<:AbstractGraph{S}}
    tensor::T
    system::L
end
