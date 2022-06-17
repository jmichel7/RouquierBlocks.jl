
# RouquierBlocks.jl

Translated from GAP3. © July 2015 — Maria Chlouveraki for the mathematics, Maria Chlouveraki and Jean Michel for the code.

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
