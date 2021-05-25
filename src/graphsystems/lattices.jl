# TODO: encode space group symmetries
const LatticeEdge{N} = Tuple{Int,Int,CartesianIndex{N}}

# N: dimensionality of the lattice, i.e. number of independent translation vectors
# D: dimensionality of vectors, typically equal to N
# M: Number of sites in unit cell

struct UnitCell{S, M, D, T<:Number}
    sites::NTuple{M, S}
    offsets::NTuple{M, SVector{D, T}}
end
UnitCell(sites::NTuple{M,S}, offsets::NTuple{M,NTuple{D,T}}) where {S, M, D, T<:Number} =
    UnitCell(sites, map(SVector, offsets))
Base.length(u::UnitCell) = length(typeof(u))
Base.length(::Type{<:UnitCell{<:Any,M}}) where M = M
@propagate_inbounds Base.getindex(u::UnitCell, i::Int) = u

Base.iterate(u::UnitCell, args...) = iterate(u.sites, args...)

function Base.show(io::IO, u::UnitCell)
    print(io, "UnitCell(")
    print(io, u.sites)
    print(io, ", (")
    for i = 1:length(u.offsets)
        i > 1 && print(io, ", ")
        print(io, u.offsets[i].data)
    end
    print(io, "))")
end

struct Lattice{S, N, M, D, U<:UnitCell{S,M},
                    B<:SMatrix{D,N}, G<:Geometry{N}} <: AbstractGraph{S}
    unitcell::U
    translationvectors::B
    bulkedges::Vector{LatticeEdge{N}}
    geometry::G
    name::String
end
Base.IteratorSize(l::Lattice) = IteratorSize(typeof(l))
Base.IteratorSize(::Type{Lattice{S,N,M,D,U,B,G}}) where {S,N,M,D,U,B,G} = IteratorSize(G)

# specific lattice constructors
# one-dimensional and quasi-one-dimensional
function chain(site, a::Number = 1; geometry::Geometry{1} = Line(), name = "chain")
    unitcell = UnitCell((site,), (SVector(zero(a)),))
    translationvectors = SMatrix{1,1}(a)
    bulkedges=[(1,1,CartesianIndex(1))]
    return Lattice(unitcell, translationvectors, bulkedges, geometry, name)
end

ladder(site, ax::Number = 1, ay::Number = ax;
        geometry::Geometry{1} = Line(), name = "ladder") =
    ladder((site,site), ax, ay; geometry = geometry, name = name)
function ladder(sites::NTuple{2}, ax::Number = 1, ay::Number = ax;
                geometry::Geometry{1} = Line(), name = "ladder")
    unitcell = UnitCell(sites, (SVector(zero(ax), zero(ay)), SVector(zero(ax), ay)))
    translationvectors = SMatrix{2,1}(ax, zero(ay))
    bulkedges=[(1,2,CartesianIndex(0)), (1,1,CartesianIndex(1)), (2,2,CartesianIndex(1))]
    return Lattice(unitcell, translationvectors, bulkedges, geometry, name)
end

# two dimensional
function square(site, a::Number = 1; geometry::Geometry{2} = Plane(), name = "square")
    unitcell = UnitCell((site,), (SVector(zero(a), zero(a)),))
    translationvectors = SMatrix{2,2}(a, zero(a), zero(a), a)
    bulkedges=[(1,1,CartesianIndex(1,0)), (1,1,CartesianIndex(0,1))]
    return Lattice(unitcell, translationvectors, bulkedges, geometry, name)
end

function rectangular(site, ax::Number, ay::Number;
        geometry::Geometry{2} = Plane(), name = "rectangular")
    unitcell = UnitCell((site,), (SVector(zero(a), zero(a)),))
    translationvectors = SMatrix{2,2}(ax, zero(ay), zero(ax), ay)
    bulkedges=[(1,1,CartesianIndex(1,0)), (1,1,CartesianIndex(0,1))]
    return Lattice(unitcell, translationvectors, bulkedges, geometry, name)
end

function oblique(site, a::Number, b::Number, θ::Number;
        geometry::Geometry{2} = Plane(), name = "oblique")
    unitcell = UnitCell((site,), (SVector(zero(a), zero(a)),))
    translationvectors = SMatrix{2,2}(a, zero(b), b*cos(θ), b*sin(θ))
    bulkedges=[(1,1,CartesianIndex(1,0)), (1,1,CartesianIndex(0,1))]
    return Lattice(unitcell, translationvectors, bulkedges, geometry, name)
end

function triangular(site, a::Number = 1; geometry::Geometry{2} = Plane(), name = "triangular")
    unitcell = UnitCell((site,), (SVector(zero(a), zero(a)),))
    sqrt3 = sqrt(oftype(a, 3))
    translationvectors = SMatrix{2,2}(a, zero(a), a/2, a*sqrt3/2)
    bulkedges=[(1,1,CartesianIndex(1,0)), (1,1,CartesianIndex(0,1)),
            (1,1,CartesianIndex(1,-1))]
    return Lattice(unitcell, translationvectors, bulkedges, geometry, name)
end

honeycomb(site, a::Number = 1; geometry::Geometry{2} = Plane(), name = "honeycomb") =
    honeycomb((site,site), a; geometry = geometry, name = name)
function honeycomb(sites::NTuple{2}, a::Number = 1;
                    geometry::Geometry{2} = Plane(), name = "honeycomb")
    unitcell = UnitCell(sites, (SVector(zero(a), zero(a)), SVector(a, zero(a))))
    translationvectors = SMatrix{2,2}(3a/2, sqrt(3)*a/2, 3a/2, -sqrt(3)a/2)
    bulkedges=[(1,2,CartesianIndex(0,0)),
                (2,1,CartesianIndex(1,0)),
                (2,1,CartesianIndex(1,0))]
    return Lattice(unitcell, translationvectors, bulkedges, geometry, name)
end

kagome(site, a::Number = 1; geometry::Geometry{2} = Plane(), name = "kagome") =
    kagome((site,site,site), a; geometry = geometry, name = name)
function kagome(sites::NTuple{3}, a::Number = 1;
        geometry::Geometry{2} = Plane(), name = "kagome")
    sqrt3 = sqrt(oftype(a, 3))
    offset1 = SVector(zero(sqrt3), zero(sqrt3))
    offset2 = SVector(a, zero(sqrt3))
    offset3 = SVector(a/2, sqrt3*a/2)
    unitcell = UnitCell(sites, (offset1, offset2, offset3))
    translationvectors = SMatrix{2,2}(2*a, zero(a), a, sqrt3*a)
    bulkedges=[(1,2,CartesianIndex(0,0)),
                (1,3,CartesianIndex(0,0)),
                (2,3,CartesianIndex(0,0)),
                (2,1,CartesianIndex(1,0)),
                (3,1,CartesianIndex(0,1)),
                (3,2,CartesianIndex(-1,1))]
    return Lattice(unitcell, translationvectors, bulkedges, geometry, name)
end

# some arbitrary normalization using integer linear combinations
function normalize_translationvectors(latticevectors::SMatrix{D,N}) where {D,N}
    vecs = collect(latticevectors[:,i] for i = 1:N)
    while true
        converged = true
        for i = 1:N
            for j = 1:N
                i == j && continue
                s = sign(dot(vecs[j],vecs[i]))
                while norm(vecs[i] - s*vecs[j]) < norm(vecs[i])
                    vecs[i] -= s*vecs[j]
                    converged = false
                end
            end
        end
        converged && break
    end
    for i = N:-1:1
        e_i = SVector{N}(ntuple(==(i), Val(N)))
        for j = 1:N
            if dot(e_i,vecs[j]) < 0
                vecs[j] = -vecs[j]
            end
        end
        sort!(vecs, by=v->dot(e_i,v), rev = true)
    end
    return SMatrix{D,N}(reshape(reinterpret(eltype(latticevectors), vecs), (D,N)))
end

# shift vectors to Wigner Seitz unit cell
function shift_to_wignerseitz(vector::SVector{D}, latticevectors::SMatrix{D,N}) where {D,N}
    R = CartesianIndices(ntuple(n->-1:1, Val(N)))
    while true
        prevvector = vector
        for I in R
            I == zero(I) && continue
            b = sum(k->I[k]*latticevectors[:,k], 1:N)
            n = round(Int, dot(b,vector)/dot(b,b), RoundNearest)
            if n > 0
                vector -= n*b
            end
        end
        prevvector == vector && break # converged
    end
    return vector
end

function _edgedirection(e::LatticeEdge{N}, l::Lattice{S,N}) where {S,N}
    i, j, I = e
    v = l.unitcell.offsets[j] - l.unitcell.offsets[i]
    for k = 1:N
        v += I[k] * l.translationvectors[:,k]
    end
    return v
end
function _edgelength(e::LatticeEdge{N}, l::Lattice{S,N}) where {S,N}
    return norm(_edgedirection(e, l))
end
_edgereverse(e::LatticeEdge) = (e[2], e[1], -e[3])


function bulkedges(l::Lattice{S,N}; order::Int = 1) where {S,N}
    if order == 0
        return LatticeEdge{N}[(i,i, zero(CartesianIndex{N})) for i = 1:length(l.unitcell)]
    else
        R = CartesianIndices(ntuple(n->(-order:order), Val(N)))
        edges = Vector{LatticeEdge{N}}(undef, 0)
        M = length(l.unitcell)
        for I in R
            k = findfirst(!iszero, Tuple(I))
            k !== nothing && I[k] < 0 && continue
            for j in Base.OneTo(M), i in Base.OneTo(I==zero(I) ? j-1 : M)
                e = (i,j,I)
                push!(edges, e)
            end
        end
        lengths = _edgelength.(edges, (l,))
        s = minimum(lengths)
        indices = findall(x->≈(x,s), lengths)
        while order > 1
            deleteat!(lengths, indices)
            deleteat!(edges, indices)
            order -= 1
            s = minimum(lengths)
            indices = findall(x->≈(x,s), lengths)
        end
        return edges[indices]
    end
end
