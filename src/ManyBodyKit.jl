module ManyBodyKit

using TensorKit
using StaticArrays
# using LinearAlgebra
using HalfIntegers
using Base.Iterators
using Base: @propagate_inbounds,
            IteratorSize, IteratorEltype, HasLength, IsInfinite, SizeUnknown, HasEltype

export Geometry, Line, HalfLine, Circle, Plane, Strip, Cylinder, Torus, ×
export UnitCell, Lattice, Site, Edge, sites, edges, position, distance, direction
export chain, ladder, square, rectangular, oblique, triangular, honeycomb, kagome

export ElementaryQuantumSystem, LocalQuantumOperator

export Qubits, Bosons, Spins

# export Lattice, Site, Edge, UnitCell, Geometry
# export linear, square, rectangular, centered_rectangular, triangular, honycomb, kagome, oblique
# export chain, ring, plane, torus, cylinder
# export position, distance, neighbors, edges, sites, bravaisindices
# export simplevec
include("graphsystems/geometry.jl")
include("graphsystems/graphs.jl")
include("graphsystems/lattices.jl")
include("graphsystems/subgraphs.jl")
include("graphsystems/iterators.jl")

include("quantumsystems/generic.jl")
include("quantumsystems/qubits.jl")
include("quantumsystems/spins.jl")
include("quantumsystems/bosons.jl")
#

# include("lattice.jl")
# include("operators.jl")
end # module
