# Leech-lattice

SageMath code for the Leech lattice $\Lambda$ and its automorphisms; in particular, matrix representatives for the 167 conjugacy classes in the Conway group $\mathrm{Co}_0 = \mathrm{Aut}(\Lambda)$.

The representatives were generated with the help of [this MathOverflow question and answer](https://mathoverflow.net/questions/338095/where-or-how-can-i-find-matrix-representatives-of-the-conjugacy-classes-of-conwa) by user Gro-Tsen.

# How to use

After downloading the two files:
start SageMath (using Python 3), load the Python script:
> load('leech_lattice_automorphisms.py')

and load the matrix representatives:
> X = load('conjugacy_class_representatives.sobj')

The Leech lattice is viewed as $\mathbb{Z}^{24}$ with the $24 \times 24$ Gram matrix $S$ determined by Curtis' Miracle Octad Generator (MOG). The Gram matrix can be called with
>Leech()

The object $X$ loaded above is a list of *LeechLatticeAutomorphisms*. Each of these is represented by a $24\times 24$ matrix $g$ satisfying $g^T S g = S$. The following methods are useful:
> g.cycle_shape()

outputs the cycle shape of $g \in X$ as a string: this is the formal product $\prod_{k \ge 1} k^{a_k}$, where $$\mathrm{det}(1 - gX) = \prod_{k \ge 1} (1 - X^k)^{a_k}.$$
> g.level()

outputs the level of $g \in X$: this is the level of the associated eta product $\eta_g(\tau) = \prod_{k \ge 1} \eta(k \tau)^{a_k}.$
> g.order()

outputs the order of $g \in X$: this is the minimal $n \in \mathbb{N}$ for which $g^n$ is the identity.

> g.fixed_lattice()
outputs the fixed lattice $\Lambda^g = \{x \in \Lambda: \, gx = x\}$ as a *LeechSubLattice* instance. This has several methods including *gram_matrix()*; *basis_matrix()*; and *orthogonal_projection(v)* which computes the orthogonal projection of a vector $v \in \Lambda$ to $\Lambda^g$. For more details please read the source code.
