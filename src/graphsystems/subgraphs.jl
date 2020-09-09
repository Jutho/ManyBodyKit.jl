struct Site{S, N, L<:Lattice{S,N}} <: Subgraph{S}
    cellindex::Int
    index::CartesianIndex{N}
    lattice::L
end
sites(s::Site) = (s,)
Base.length(s::Site) = 1
Base.getindex(s::Site) = s.lattice.unitcell.sites[s.cellindex]

struct Edge{S, N, L<:Lattice{S,N}} <: Subgraph{S}
    site1::Site{S,N,L}
    site2::Site{S,N,L}
end
sites(e::Edge) = (e.site1, e.site2)
Base.length(e::Edge) = 2
Base.reverse(e::Edge) = Edge(e.site2, e.site1)
Base.iterate(e::Edge, args...) = iterate((e.site1, e.site2), args...)
Base.getindex(e::Edge, i) = i == 1 ? e.site1 : (i == 2 ? e.site2 : throw(BoundsError(e, i)))

function position(s::Site)
    l = s.lattice
    u = l.unitcell
    g = l.geometry
    c = u.offsets[s.cellindex]
    for i = 1:length(s.index)
        c += s.index[i] * l.translationvectors[:,i]
    end
    return c
end
position(e::Edge) = (position(e.site1) + position(e.site2))//2
direction(e::Edge) = (position(e.site2) - position(e.site1))

function distance(s1::Site, s2::Site)
    @assert s1.lattice === s2.lattice
    return norm(position(s1) - position(s2))
end
function neighbors(s::Site; order::Int = 1)
    order == 0 && return [s]
    @assert order > 0
    lattice = s.lattice
    geom = lattice.geometry
    neighbors = Vector{typeof(s)}(undef, 0)
    I = s.index
    for (i,j,ΔI) in bulkedges(lattice; order = order)
        J = I + ΔI
        if i == s.cellindex && J ∈ lattice.geometry
            push!(neighbors, Site(j, J, lattice))
        end
        J = I - ΔI
        if j == s.cellindex && J ∈ lattice.geometry
            push!(neighbors, Site(i, J, lattice))
        end
    end
    return neighbors
end
edges(s::Site; order = 1) = Edge.((s,), neighbors(s; order = order))
