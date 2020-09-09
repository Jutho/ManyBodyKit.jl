abstract type Extent end

struct Infinite <: Extent
end

struct HalfInfinite <: Extent
    first::Int
    step::Int
    function HalfInfinite(first = 0, step = 1)
        abs(step) == 1 || throw(ArgumentError("step of HalfInfinite can only be +1 or -1"))
        new(first,step)
    end
end

struct Open <: Extent
    range::UnitRange{Int}
end
Open(L::Int) = Open(1:L)
Base.length(e::Open) = length(e.range)

struct Closed <: Extent
    length::Int
end
Base.length(e::Closed) = e.length

Base.in(i::Int, e::Open) = i in e.range
Base.in(i::Int, e::Closed) = true
Base.in(i::Int, e::HalfInfinite) = (i - e.first)*e.step >= 0
Base.in(i::Int, e::Infinite) = true

Base.mod(i::Int, e::Open) = i in e ? i : throw(BoundsError(e, i))
Base.mod(i::Int, e::Closed) = mod(i, e.length)
Base.mod(i::Int, e::HalfInfinite) = i in e ? i : throw(BoundsError(i, e))
Base.mod(i::Int, e::Infinite) = i

Base.iterate(e::Open, args...) = iterate(e.range, args...)
Base.iterate(e::Closed, i = 0) = i == e.length ? nothing : (i, i+1)
Base.iterate(e::HalfInfinite, i = e.first) = (i, i+e.step)
Base.iterate(e::Infinite, i = 0) = (i, i <= 0 ? -i+1 : -i)

Base.IteratorSize(::Type{<:Union{Open,Closed}}) = Base.HasLength()
Base.IteratorSize(::Type{<:Union{Infinite,HalfInfinite}}) = Base.IsInfinite()

Base.IteratorEltype(::Type{<:Extent}) = Base.HasEltype()
Base.eltype(::Type{<:Extent}) = Int

const Twist{N} = NTuple{N,CartesianIndex{N}}

struct Geometry{N, G<:NTuple{N,Extent}}
    extents::G
    twists::Twist{N}
    function Geometry{N, G}(extents::G, twists::Twist{N}) where {N, G}
        for i = 1:N
            if extents[i] isa Closed
                for k = 1:i
                    twists[i][k] == 0 ||
                        throw(ArgumentError("twist in direction $i should have zero as its first $(i) components."))
                end
            else
                twists[i] == zero(twists[i]) ||
                    throw(ArgumentError("twist should be zero for non-`Closed` directions."))
            end
        end
        # TODO: better check twists? more flexible twists?
        return new(extents, twists)
    end
end
Geometry(extents::NTuple{N,Extent}) where N =
    Geometry(extents, ntuple(n->zero(CartesianIndex{N}), Val(N)))
Geometry(extents::NTuple{N,Extent}, twists::Twist{N}) where N =
    Geometry{N,typeof(extents)}(extents, twists)

function Base.iterate(g::Geometry, args...)
    next = iterate(Base.Iterators.product(g.extents...), args...)
    next === nothing && return next
    return CartesianIndex(next[1]), next[2]
end
Base.IteratorSize(::Type{Geometry{N,G}}) where {N,G} =
    Base.IteratorSize(Iterators.ProductIterator{G})
Base.size(g::Geometry) = IteratorSize(g) != IsInfinite() ? map(length, g.extents) :
    throw(ArgumentError("no size for infinite geometry"))
Base.length(g::Geometry) = IteratorSize(g) != IsInfinite() ? prod(length, g.extents) :
    throw(ArgumentError("no length for infinite geometry"))
Base.IteratorEltype(::Type{<:Geometry}) = Base.HasEltype()
Base.eltype(::Type{<:Geometry{N}}) where N = CartesianIndex{N}
Base.eltype(g::Geometry) = eltype(typeof(g))

Base.in(I::CartesianIndex{N}, g::Geometry{N}) where N = all(in.(Tuple(I), g.extents))
function Base.mod(I::CartesianIndex{N}, g::Geometry{N}) where N
    I in g || throw(BoundsError(g, I))
    for n = 1:N
        e = g.extents[n]
        if e isa Closed
            T_n = CartesianIndex(ntuple(k->ifelse(n==k, length(e), 0), Val(N)))
            while I[n] >= length(e)
                I -= T_n
                I += g.twists[n]
            end
            while I[n] < 0
                I += T_n
                I -= g.twists[n]
            end
        end
    end
    return I
end

# specific geometries
const HalfLine = Geometry{1,Tuple{HalfInfinite}}
const Line = Geometry{1,<:Tuple{Union{Infinite,Open}}}
const Circle = Geometry{1,Tuple{Closed}}
const Plane = Geometry{2,<:NTuple{2,Union{Infinite,Open}}}
const Strip = Geometry{2,Tuple{Open,Infinite}}
const Cylinder = Geometry{2,Tuple{Closed,Union{Infinite,Open}}}
const Torus = Geometry{2,Tuple{Closed,Closed}}

HalfLine() = Geometry((HalfInfinite(),))
Line() = Geometry((Infinite(),))
Line(L::Integer) = Geometry((Open(L),))
Circle(circumference::Integer) = Geometry((Closed(circumference),))
Plane(width::Integer, height::Integer) = Geometry((Open(width), Open(height)))
Plane() = Geometry((Infinite(), Infinite()))
Strip(width::Integer) = Geometry((Open(width), Infinite()))
Cylinder(circumference::Integer, length::Integer) =
    Geometry((Closed(circumference), Open(length)))
Cylinder(circumference::Integer; twist::Integer = 0) =
    Geometry((Closed(circumference::Integer), Infinite()),
                CartesianIndex.(((0,twist), (0,0))))
Torus(width::Integer, height::Integer; twist::Int = 0) =
    Geometry((Closed(width), Closed(height)), CartesianIndex.(((0,twist), (0,0))))

# combine geometries
function TensorKit.:×(g1::Geometry{N1,G1}, g2::Geometry{N2,G2}) where {N1,N2,G1,G2}
    extents = (g1.extents..., g2.extents...)
    zeroN1 = zero(CartesianIndex{N1})
    zeroN2 = zero(CartesianIndex{N2})
    twists = ntuple(Val(N1+N2)) do i
        if i <= N1
            return CartesianIndex(Tuple(g1.twists[i])..., Tuple(zeroN2)...)
        else
            return CartesianIndex(Tuple(zeroN1)..., Tuple(g2.twists[i-N1])...)
        end
    end
    return Geometry(extents, twists)
end

function Base.show(io::IO, g::Geometry{N}) where N
    if N == 1
        e, = g.extents
        if e isa Infinite
            return Base.print(io, "Line()")
        elseif e isa HalfInfinite
            return Base.print(io, "HalfLine($(e.first), $(e.step))")
        elseif e isa Open
            return Base.print(io, "Line($(e.range))")
        elseif e isa Closed
            return Base.print(io, "Circle($(e.length))")
        end
    elseif N == 2
        e1, e2 = g.extents
        if e2 isa Infinite
            if e1 isa Infinite
                return print(io, "Plane()")
            elseif e1 isa Open
                return print(io, "Strip($(e1.range))")
            elseif e1 isa Closed
                print(io, "Cylinder(")
                print(io, e1.length)
                if !all(iszero, g.twists)
                    t = g.twists[1][2]
                    print(io, "; twist = $t")
                end
                return print(io, ")")
            end
        elseif e2 isa Open && e1 isa Open
            return print(io, "Plane($(e1.range), $(e2.range))")
        elseif e2 isa Closed && e1 isa Closed
            print(io, "Torus($(e1.length), $(e2.length)")
            if !all(iszero, g.twists)
                t = g.twists[1][2]
                print(io, "; twist = $t")
            end
            return print(io, ")")
        end
    end
    Base.print(io, "Geometry(")
    Base.print(io, g.extents)
    if !all(iszero, g.twists)
        Base.print(io, ", ")
        Base.print(io, g.twists)
    end
    Base.print(")")
end

#
# function Base.show(io::IO, g::Geometry{N}) where N
#     if N == 1
#         e, = g.extents
#         if e isa Infinite
#             Base.print(io, "line()")
#         elseif e isa HalfInfinite
#             Base.print(io, "halfline($(e.first), $(e.step))")
#         elseif e isa Open
#             Base.print(io, "line($(e.range))")
#         elseif e isa Closed
#             Base.print(io, "circle($(e.length))")
#         end
#     elseif N == 2
#         e1, e2 = g.extents
#     end
# end
