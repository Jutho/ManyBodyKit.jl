using TensorKit, ManyBodyKit


### Classical things
## uniform Ising
using ManyBodyKit.Classical
physical = classical(2) # 2 local dof?
lattice = square(physical)

ham = J * ∑ σᵢσⱼ for (i,j) in neighbors(lattice)


## non-uniform Ising (v1)
physical = classical(2) 
lattice = square(physical)

vertical = L ∑ σᵢσⱼ for (i,j) in neighbors(lattice, dir=1)
horizontal = K ∑ σᵢσⱼ for (i,j) in neighbors(lattice, dir=2)

ham = vertical + horizontal


## non-uniform Ising (v2)
physical = classical(2)
lattice = rectangular(physical, ratio) # non-equivalent spacing in 2 dimensions

J(i,j) = cst / dist(i,j)
ham = ∑ J(i,j) * σᵢσⱼ for (i,j) in neighbors(lattice)


### Kagome things
## option 1: Kagome lattice knows how to determine neighbours, next-to-nearest, etc?

## option 2: Kagome lattice <: (oblique) Bravais lattice containing 2 sites?
