using TensorKit
using ManyBodyKit
using ManyBodyKit.Bosons

physical = boson(4);
lattice = chain(physical);

hop1 = ∑(b⁺(i)⊗b⁻(j) for (i,j) in neighbors(lattice));
hop2 = ∑(b⁻(i)⊗b⁺(i+1) for i in sites(lattice));
ham = ∑(hop1,hop2,(n²(i) for i in sites(lattice)));
