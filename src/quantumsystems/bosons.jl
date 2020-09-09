module Bosons
    using TensorKit, LinearAlgebra, ..ManyBodyKit

    export boson, annihilator, creator, density, density2, hopping
    export b⁻, b⁺, n, n²

    # spaces
    boson(n::Int, ::Type{Trivial} = Trivial) = ℂ^(n+1)
    boson(n::Int, ::Type{ℤ₂}) = Z₂Space(0=>length(0:2:n), 1=>length(1:2:n))
    boson(n::Int, ::Type{U₁}) = U₁Space(U₁(k)=>1 for k=0:n)

    # operators
    for (op, s) in ((:annihilator, :b⁻), (:creator, :b⁺), (:density, :n), (:density2, :n²))
        @eval const $op = LocalQuantumOperator{$(QuoteNode(s))}()
        @eval const $s = $op
    end
    const hopping = LocalQuantumOperator{:bosonhopping}

    for op in (:annihilator, :creator, :density, :hopping)
        @eval $op(V::ElementarySpace) = $op(Float64, V)
    end

    function annihilator(T::Type{<:Number}, V::ComplexSpace)
        data = zeros(T, dim(V), dim(V))
        for n = 1:dim(V)-1
            data[n,n+1] = sqrt(n)
        end
        return TensorMap(data, V, V)
    end

    function creator(T::Type{<:Number}, V::ComplexSpace)
        data = zeros(T, dim(V), dim(V))
        for n = 1:dim(V)-1
            data[n+1,n] = sqrt(n)
        end
        return TensorMap(data, V, V)
    end

    function density(T::Type{<:Number}, V::ComplexSpace)
        data = zeros(T, dim(V), dim(V))
        for n = 1:dim(V)-1
            data[n+1,n+1] = n
        end
        return TensorMap(data, V, V)
    end

    function annihilator(T::Type{<:Number}, V::U₁Space)
        t = TensorMap(zeros, T, V ⊗ U₁Space(1=>1), V)
        for (n, b) in blocks(t)
            size(b) == (1,1) || throw(DomainError("not a valid boson space"))
            b[1,1] = sqrt(n.charge)
        end
        return t
    end
    creator(T, V::U₁Space) = annihilator(T, V)'

    function density(T, V::U₁Space)
        t = TensorMap(zeros, T, V ⊗ U₁Space(1=>1), V)
        for (n, b) in blocks(t)
            size(b) == (1,1) || throw(DomainError("not a valid boson space"))
            isinteger(n.charge) || throw(DomainError("half-integer boson charge impossible"))
            b[1,1] = n.charge
        end
        return t
    end

    # two-site operators
    hopping(T, V::ProductSpace{ComplexSpace,2}) =
        creator(T, V[1]) ⊗ annihilator(T, V[2]) + annihilator(T, V[1]) ⊗ creator(T, V[2])

    function hopping(T, V::ProductSpace{U₁Space,2})
        @tensor h[a b; a' b'] :=
            creator(T, V[1])[a,a',s] * annihilator(T, V[2])[b,s,b'] +
            annihilator(T, V[1])[a,s,a'] * creator(T, V[2])[b,b',s]
        return h
    end
end
