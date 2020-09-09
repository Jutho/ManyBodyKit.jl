sites(l::Lattice) = SiteIterator(l)
edges(l::Lattice; order::Int = 1) = EdgeIterator(l, order)

struct SiteIterator{L<:Lattice}
    lattice::L
end
struct EdgeIterator{L<:Lattice}
    lattice::L
    order::Int
end
Base.IteratorSize(::SiteIterator{L}) where L = Base.IteratorSize(L)
Base.IteratorSize(::EdgeIterator{L}) where L = Base.IteratorSize(L)

Base.size(iter::SiteIterator) =
    (length(iter.lattice.unitcell), size(iter.lattice.geometry)...)

function Base.iterate(iter::SiteIterator, args...)
    lattice = iter.lattice
    geometry = lattice.geometry
    unitcell = lattice.unitcell
    r = Base.OneTo(length(unitcell))
    next = iterate(Base.Iterators.product(r, geometry), args...)
    next === nothing && return nothing
    (i, I), state = next
    return Site(i, I, lattice), state
end

#
# (Site(site,i,I,l) for ((i,site), I) in product(l.unitcell, l.geometry))
