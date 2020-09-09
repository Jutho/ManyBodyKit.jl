abstract type AbstractGraph{S} end
abstract type Subgraph{S} <: AbstractGraph{S} end

# Subgraph types should have a finite number of sites
Base.IteratorSize(::Type{<:Subgraph}) = Base.HasLength()
