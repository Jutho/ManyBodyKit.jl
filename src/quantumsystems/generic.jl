const ElementaryQuantumSystem = EuclideanSpace{ℂ}
const CompositeQuantumSystem{S} = CompositeSpace{S} where {S<:ElementaryQuantumSystem}
const QuantumSystem{S} = Union{S,CompositeSpace{S}} where {S<:ElementaryQuantumSystem}

const QuantumLatticeSystem{S} = Lattice{S} where {S<:ElementaryQuantumSystem}
const QuantumSubgraph{S} = Subgraph{S} where {S<:ElementaryQuantumSystem}
const QuantumSite{S,N,L} = Site{S,N,L} where {S<:ElementaryQuantumSystem,N,L<:Lattice{S,N}}

TensorKit.space(s::QuantumSite) = s[]
TensorKit.space(g::QuantumSubgraph) =
    IteratorSize(sites(g)) !== IsInfinite() ? ⊗(map(space, sites(g))...) :
        throw(DomainError(g, "The space of an infinite graph cannot be computed"))

abstract type QuantumOperator end
struct LocalQuantumOperator{s} <: QuantumOperator end

(O::LocalQuantumOperator)(g::Subgraph) =
    MultiSiteOperator(O(space(g)), tuple(sites(g)...))

struct MultiSiteOperator{K,S,T<:TensorMap{S,K,K},N,L<:AbstractGraph{S}} <: QuantumOperator
    tensor::T
    sites::NTuple{K,Site{S,N,L}}
end

Base.:*(a::Number, o::MultiSiteOperator) = MultiSiteOperator(a*o.tensor, o.sites)

# operators
export id
id_symbol = Symbol("𝟙")
const id = LocalQuantumOperator{id_symbol}()
id(s::ElementarySpace) = id(Bool, s)
id(T::Type{<:Number}, s::ElementarySpace) = TensorMap(UniformScaling(one(T)), s, s)
@static if VERSION > v"1.3-"
    @eval const $id_symbol = id # 𝟙 does not work as identifier in earlier Julia
end

# function generators(T::Type{<:Number}, V::U₁Space)
#     t = TensorMap(zeros, T, V, V ⊗ U₁Space(0=>1))
#     for (n,b) in blocks(t)
#         copyto!(b, n.charge * I)
#     end
#     return t
# end
# function generators(T::Type{<:Number}, V::SU₂Space)
#     W = SU₂Space(1=>1)
#     t = TensorMap(ones, T, V, V ⊗ W)
#     for (s, b) in blocks(t)
#         copyto!(b, sqrt(s.j*(s.j+1))*I)
#     end
#     return t
# end
