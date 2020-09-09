module Qubits
    using TensorKit, LinearAlgebra, ..ManyBodyKit

    export qubit, qutrid, qudit, sigmax, sigmay, sigmaz, sigmam, sigmap, paulis
    export σˣ, σʸ, σᶻ, σ⁺, σ⁻

    # spaces
    qudit(d::Int) = ℂ^d
    qubit() = qudit(2)
    qutrit() = qudit(3)

    for (op, s) in ((:sigmax, :σˣ), (:sigmay, :σʸ), (:sigmaz, :σᶻ), (:sigmap, :σ⁺), (:sigmam, :σ⁻),)
        @eval const $op = LocalQuantumOperator{$(QuoteNode(s))}()
        @eval const $s = $op
    end

    # operators
    for op in (:sigmax, :sigmay, :sigmaz, :sigmap, :sigmam, :paulis)
        @eval $op() = $op(ℂ^2)
        @eval $op(T::Type{<:Number}) = $op(T, ℂ^2)
        T = op == :sigmay ? Complex{Int} : Int
        @eval $op(V::ElementarySpace) = $op($T, V)
    end

    function sigmax(T::Type{<:Number}, V::ComplexSpace)
        dim(V) == 2 || throw(DomainError(V, "Pauli's σˣ only exists on ComplexSpace(2)"))
        data = T[0 1; 1 0]
        return TensorMap(data, V, V)
    end
    function sigmay(T::Type{<:Number}, V::ComplexSpace)
        dim(V) == 2 || throw(DomainError(V, "Pauli's σʸ only exists on ComplexSpace(2)"))
        data = T[0 -im; im 0]
        return TensorMap(data, V, V)
    end
    function sigmaz(T::Type{<:Number}, V::ComplexSpace)
        dim(V) == 2 || throw(DomainError(V, "Pauli's σᶻ only exists on ComplexSpace(2)"))
        data = T[1 0; 0 -1]
        return TensorMap(data, V, V)
    end
    function sigmap(T::Type{<:Number}, V::ComplexSpace)
        dim(V) == 2 || throw(DomainError(V, "Pauli's σᶻ only exists on ComplexSpace(2)"))
        data = T[0 1; 0 0]
        return TensorMap(data, V, V)
    end
    function sigmam(T::Type{<:Number}, V::ComplexSpace)
        dim(V) == 2 || throw(DomainError(V, "Pauli's σᶻ only exists on ComplexSpace(2)"))
        data = T[0 0; 1 0]
        return TensorMap(data, V, V)
    end

    function paulis(T::Type{<:Number}, V::ComplexSpace)
        dim(V) == 2 ||
            throw(DomainError(V, "Pauli operators only exists on ComplexSpace(2)"))
        data = hcat(T[0 1; 1 0], T[0 -im; im 0], T[1 0; 0 -1])
        return TensorMap(data, V, V ⊗ ℂ^3)
    end
    #
    # const σˣ = sigmax
    # const σʸ = sigmay
    # const σᶻ = sigmaz
    # const σ⁺ = sigmap
    # const σ⁻ = sigmam
end
