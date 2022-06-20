
# RouquierBlocks.jl

Computation of *Rouquier blocks* and *essential hyperplanes* of cyclotomic Hecke algebras of complex reflection groups. Translated from GAP3. © July 2015 — [Maria Chlouveraki](https://chlouveraki.perso.math.cnrs.fr) for the mathematics, Maria Chlouveraki and [Jean Michel](https://webusers.imj-prg.fr/~jean.michel/) for the code.

## Installation

```julia
julia> using Pkg

julia> Pkg.add(url="https://github.com/jmichel7/RouquierBlocks.jl")
```

This packages requires [Gapjm.jl](https://github.com/jmichel7/Gapjm.jl).

## Usage

```julia
julia> using RouquierBlocks, Gapjm

# Complex reflection group G4
julia> W = complex_reflection_group(4)
G₄

# Essential hyperperplanes and Rouquier blocks
julia> rouquier_blocks(W ; names=true, limit=true)
7-element Vector{Vector{Vector}}:
 [[0, 0, 0], [["φ₁‚₀"], ["φ₁‚₄"], ["φ₁‚₈"], ["φ₂‚₅"], ["φ₂‚₃"], ["φ₂‚₁"], ["φ₃‚₂"]]]
 [[0, 1, -1], [["φ₁‚₀"], ["φ₁‚₄", "φ₁‚₈", "φ₂‚₅"], ["φ₂‚₃", "φ₂‚₁"], ["φ₃‚₂"]]]
 [[1, -2, 1], [["φ₁‚₀"], ["φ₁‚₄", "φ₂‚₃", "φ₃‚₂"], ["φ₁‚₈"], ["φ₂‚₅"], ["φ₂‚₁"]]]
 [[1, -1, 0], [["φ₁‚₀", "φ₁‚₄", "φ₂‚₁"], ["φ₁‚₈"], ["φ₂‚₅", "φ₂‚₃"], ["φ₃‚₂"]]]
 [[1, 0, -1], [["φ₁‚₀", "φ₁‚₈", "φ₂‚₃"], ["φ₁‚₄"], ["φ₂‚₅", "φ₂‚₁"], ["φ₃‚₂"]]]
 [[1, 1, -2], [["φ₁‚₀"], ["φ₁‚₄"], ["φ₁‚₈", "φ₂‚₁", "φ₃‚₂"], ["φ₂‚₅"], ["φ₂‚₃"]]]
 [[2, -1, -1], [["φ₁‚₀", "φ₂‚₅", "φ₃‚₂"], ["φ₁‚₄"], ["φ₁‚₈"], ["φ₂‚₃"], ["φ₂‚₁"]]]

# The entries in this list are pairs of a list giving the coefficients of an essential
# hyperplane and a list giving the Rouquier blocks on this hyperplane. The basis of the
# parameter space is indexed by (Ω,j) with 1 ≤ j ≤ e_Ω and the Ω are in the same order as 
# in hyperplane_orbits(W).
# The first entry [0,0,...,0] is for the zero-"hyperplane", which means generic parameters.
# The optional arguments are for printing characternames: with names=false (default), the 
# indices as in the character table of W are returned.

# From the above data we can obtain the rouquier blocks at arbitrary parameters.
# We create a 1-cyclotomic specializion of the Hecke algebra for the parameter [1,0,-1], 
# which lies (generically) on the essential hyperplane [1,-2,1]. To this end, we need to
# specify a variable x.
julia> @Mvp x
Mvp{Int64}: x

# 1-cyclotomic Hecke algebra at [1,0,-1]
julia> H = hecke(W,[ [E(3)^0 * x^1, E(3)^1 * x^0, E(3)^2 * x^(-1)] ])
hecke(G₄,Vector{Mvp{Cyc{Int64}, Int64}}[[x, ζ₃, ζ₃²x⁻¹]])

julia> rouquier_blocks(H ; names=true, limit=true)
5-element Vector{Vector{String}}:
 ["φ₁‚₀"]
 ["φ₁‚₄", "φ₂‚₃", "φ₃‚₂"]
 ["φ₁‚₈"]
 ["φ₂‚₅"]
 ["φ₂‚₁"]
```

## References

1. Chlouveraki, M. (2009). *Blocks and families for cyclotomic Hecke algebras*. Springer-Verlag, Berlin. https://arxiv.org/abs/0807.1476
2. Bonnafé, C. & Rouquier, R. (2017). *Cherednik algebras and Calogero-Moser cells*. [https://arxiv.org/abs/1708.09764](https://arxiv.org/abs/1708.09764)

<!--
A  1-cyclotomic Hecke  algebra for  the complex  reflection group `W` is an Hecke algebra `H` whose `j`-th parameter for the `i`-th generator (of order `e`)  of `W` is of the form  `ζₑʲ x^mᵢ,ⱼ` for some rational numbers `mᵢ,ⱼ`; thus such an algebra specializes to the group algebra for x->1.

In this module `x` must be `Mvp(:x)`.

A  tool  to  determine  the  Rouquier  blocks  of  `H`  are  the "essential hyperplanes"  which are integral  linear forms (operating  on the variables `mᵢ,ⱼ`)  determined by the Schur elements of the generic algebra associated to   `H`.  For  each  essential  hyperplane  `h`  there  is  an  associated 1-cyclotomic  algebra  `A_h`  whose  `mᵢ,ⱼ`  annihilate  `h`  and  no other essential  hyperplane. Let us call `h`-blocks the Rouquier blocks of `A_h`. Then  the Rouquier  blocks of  `H` are  the lcm  of the  `h`-blocks for `h` running  over the hyperplanes  annihilating the `mᵢ,ⱼ`  of `H`. In the case where  the `mᵢ,ⱼ`  annihilate no  hyperplane we  get the 0-blocks.

<a id='RouquierBlocks.rouquier_blocks' href='#RouquierBlocks.rouquier_blocks'>#</a>
**`RouquierBlocks.rouquier_blocks`** &mdash; *Function*.



The  Rouquier blocks  of a  1-cyclotomic algebra  H is the finest partition coarser than h-blocks for all hyperplanes h annihilated by H's parameters.

<a id='RouquierBlocks.RouquierBlockData' href='#RouquierBlocks.RouquierBlockData'>#</a>
**`RouquierBlocks.RouquierBlockData`** &mdash; *Function*.



`RouquierBlockData(W)`

returns  a list of [essential hyperplane h, corresponding h-blocks] for the complex reflection group `W`.

`h`  is represented  as a  list of  integers of  same length as the list of parameters  for  the  Hecke  algebra  of  `W`. `h`-blocks is a partition of `1:nconjugacy_classes(W)`.

The  first entry in the result list has `h=[0,...,0]` and the corresponding `h`-blocks are the `0`-blocks.
-->