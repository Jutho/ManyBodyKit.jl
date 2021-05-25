using TensorKit
using ManyBodyKit
using ManyBodyKit.Spins

len = 10;
ph1 = spin(1//2);
ph2 = spin(1);
lattice = chain((ph1,ph2),geometry = Circle(len));

prefactors = rand(len);
ham = ∑(prefactors[position(i)]*Sz(i)⊗Sz(j) for (i,j) in neighbors(lattice));
