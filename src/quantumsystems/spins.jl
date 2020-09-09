module Spins
    using TensorKit, LinearAlgebra, HalfIntegers, ..ManyBodyKit

    export spin, Sx, Sy, Sz, Sp, Sm
    export Sˣ, Sʸ, Sᶻ, S⁺, S⁻

    # spaces
    spin(J::Number, S::Type{<:Sector} = Trivial) =
        ishalfinteger(J) ? spin(HalfInteger(J), S) :
            throw(ArgumentError("spin quantum number should be half integer"))
    spin(J::HalfInteger, ::Type{Trivial} = Trivial) = ℂ^(twice(J)+1)
    spin(J::HalfInteger, ::Type{U₁}) = U₁Space(U₁(Jz)=>1 for Jz in -J:J)
    spin(J::HalfInteger, ::Type{CU₁}) = CU₁Space(CU₁(Jz)=>1 for Jz in J:-1:0)
    spin(J::HalfInteger, ::Type{SU₂}) = SU₂Space(SU₂(J)=>1)

    # operators
    for (op, s) in ((:Sx, :Sˣ), (:Sy, :Sʸ), (:Sz, :Sᶻ), (:Sp, :S⁺), (:Sm, :S⁻),)
        @eval const $op = LocalQuantumOperator{$(QuoteNode(s))}()
        @eval const $s = $op
    end
    for op in (:Sx, :Sy, :Sz, :Sp, :Sm)
        T = op == :Sy ? ComplexF64 : Float64
        @eval $op(V) = $op($T, V)
    end

    function Sx(T::Type{<:Number}, V::ComplexSpace)
        d = dim(V)
        j = (d-1)//2
        data = T[1/2*sqrt(((a==b+1) + (a==b-1))*((j+1)*(a+b-1)-a*b)) for a in 1:d, b in 1:d]
        return TensorMap(data, V, V)
    end
    function Sy(T::Type{<:Number}, V::ComplexSpace)
        d = dim(V)
        j = (d-1)//2
        data = T[im//2*sqrt(((a==b+1) - (a==b-1))*((j+1)*(a+b-1)-a*b)) for a in 1:d, b in 1:d]
        return TensorMap(data, V, V)
    end
    function Sz(T::Type{<:Number}, V::ComplexSpace)
        d = dim(V)
        j = (d-1)//2
        data = T[(a==b)*(j+1-a) for a in 1:d, b in 1:d]
        return TensorMap(data, V, V)
    end
    function Sp(T::Type{<:Number}, V::ComplexSpace)
        d = dim(V)
        j = (d-1)//2
        data = T[sqrt((a==b-1)*((j+1)*(a+b-1)-a*b)) for a in 1:d, b in 1:d]
        return TensorMap(data, V, V)
    end
    function Sm(T::Type{<:Number}, V::ComplexSpace)
        d = dim(V)
        j = (d-1)//2
        data = T[sqrt((a==b+1)*((j+1)*(a+b-1)-a*b)) for a in 1:d, b in 1:d]
        return TensorMap(data, V, V)
    end
    function Sz(T::Type{<:Number}, V::U₁Space)
        op = TensorMap(zeros, V, V)
        for (s, b) in blocks(op)
            size(b) == (1,1) || throw(DomainError("not a valid spin space"))
            b[1,1] = s.charge
        end
        return op
    end
    function Sm(T::Type{<:Number}, V::U₁Space)
        op = TensorMap(zeros, V, V * U₁Space(-1=>1))
        d = dim(V)
        j = (d-1)//2
        for (s, b) in blocks(op)
            size(b) == (1,1) || throw(DomainError("not a valid spin space"))
            m = s.charge
            b[1,1] = sqrt(j*(j+1) - m*(m+1))
        end
        return op
    end
    Sp(T::Type{<:Number}, V::U₁Space) = Sm(T, V)'

    # # two-site operators
    # for (O, o) in ((:Sxx, :Sx), (:Syy, :Sy), (:Szz, :Sz))
    #     @eval $O(T::Type{<:Number}, V::ProductSpace{ComplexSpace,2}) =
    #         $o(T, V[1]) ⊗ $o(T, V[2])
    # end
end
