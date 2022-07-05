"""
translated from gap3.
© July 2015 --- Maria Chlouveraki for the mathematics,
                Maria Chlouveraki and Jean Michel for the code

A  1-cyclotomic Hecke  algebra for  the complex  reflection group `W` is an
Hecke algebra `H` whose `j`-th parameter for the `i`-th generator (of order
`e`)  of `W` is of the form  `ζₑʲ x^mᵢ,ⱼ` for some rational numbers `mᵢ,ⱼ`;
thus such an algebra specializes to the group algebra for x->1.

In this module `x` must be `Mvp(:x)`.

A  tool  to  determine  the  Rouquier  blocks  of  `H`  are  the "essential
hyperplanes"  which are integral  linear forms (operating  on the variables
`mᵢ,ⱼ`)  determined by the Schur elements of the generic algebra associated
to   `H`.  For  each  essential  hyperplane  `h`  there  is  an  associated
1-cyclotomic  algebra  `A_h`  whose  `mᵢ,ⱼ`  annihilate  `h`  and  no other
essential  hyperplane. Let us call `h`-blocks the Rouquier blocks of `A_h`.
Then  the Rouquier  blocks of  `H` are  the lcm  of the  `h`-blocks for `h`
running  over the hyperplanes  annihilating the `mᵢ,ⱼ`  of `H`. In the case
where  the `mᵢ,ⱼ`  annihilate no  hyperplane we  get the 0-blocks.
"""
module RouquierBlocks
using Gapjm
export rouquier_blocks, rouquier_blocks, generic_schur_elements

# largest v∈ ℤ such that p^-v c is a p-algebraic integer.
function LaurentPolynomials.valuation(c::Cyc,p::Integer)
  e=conjugates(c)*big(1)
  valuation(prod(e)[0],p)//length(e)
end

function LaurentPolynomials.valuation(c::Mvp,p::Integer)
  if iszero(c) return 0 end
  minimum(valuation.(coefficients(c),p))
end

ispInt(c::Cyc,p::Integer)=mod(denominator(c),p)!=0
ispInt(c::Mvp,p::Integer)=all(x->ispInt(x,p),coefficients(c))
ispInt(c::Unknown,p::Integer)=true

"""
    function generic_schur_elements(W::ComplexReflectionGroup)

generic_schur_elements(W)  returns  the  Schur  elements for Hecke(W) with
parameters  ζₑⁱ xⱼ,ᵢ (a variant of  the generic algebra which specializes
to  the group  algebra for  xᵢ,ⱼ->1). The  result is  a named  tuple with
fields  .coeff, .mon representing the leading monomial, and a field .vcyc
which  is a list  of namedtuples with  fields .coeff, .mon representing a
monomial  m  and  a  field  .pol  holding  a  cycpol such that the record
represents pol(m)
"""
function generic_schur_elements(W::ComplexReflectionGroup)
  vars="xyz" # for irreducible groups at most 3 orbits
  hh=hyperplane_orbits(W)
  vv=[Symbol.(vars[i],0:h.order-1) for (i,h) in enumerate(hh)]
  H=hecke(W,Tuple(Mvp.(v).*E.(length(v),0:length(v)-1) for v in vv))
  vnames=vcat(vv...)
  map(FactorizedSchurElements(H))do s
    function montovec(mon)
      mon=first(monomials(mon))
      res=map(x->0//1,vnames)
      res[indexin(variables(mon),vnames)]=collect(powers(mon))
      res
    end
    res=(vars=Mvp.(vnames),vcyc=map(p->(pol=p.pol,
      mon=montovec(p.monomial),
      coeff=first(coefficients(p.monomial))),s.vcyc),
      coeff=first(coefficients(s.factor)), mon=montovec(s.factor))
    sort!(res.vcyc, by=x->[x.pol, x.mon, x.coeff])
    res
  end
end

"""
    function rouquier_blocks(W::ComplexReflectionGroup)

returns  a list of [essential hyperplane h, corresponding h-blocks] for the
complex reflection group `W`.

`h`  is represented  as a  list of  integers of  same length as the list of
parameters  for  the  Hecke  algebra  of  `W`. `h`-blocks is a partition of
`1:nconjugacy_classes(W)`.

The  first entry in the result list has `h=[0,...,0]` and the corresponding
`h`-blocks are the `0`-blocks.
"""
function rouquier_blocks(W::ComplexReflectionGroup ; names=false, namesargs...)
  bl0=nothing
  NRPARA=5 # how many random algebras A_h to consider
  sch=generic_schur_elements(W)
  hplanes=vcat(map(x->map(y->y.mon,x.vcyc),sch)...)
  hplanes=map(v->Int.(v*lcm(denominator.(v))),hplanes)
  hplanes=unique(map(v->div.(v,gcd(numerator.(v))),hplanes))
  sort!(hplanes)
  pushfirst!(hplanes,zero(hplanes[1])) # [0,..,0]+essential hplanes
  InfoChevie("#I ",length(hplanes)," hplanes\n")
  return map(hplanes)do h # for each hplane h return [h,Rouquier blocks of A_h]
    hh=filter(k->!iszero(k) && k!=h, hplanes)
    m=lnullspaceInt(hcat(h))
    para=Vector{Int}[]
    while length(para)<NRPARA
      if isempty(hh) && h==[1,-1] push!(para,[1,1])
      else
        p=rand(-2*size(m,1):2*size(m,1),size(m,1))'*m
        if !(0 in toM(hh)*p') push!(para, vec(div.(p,gcd(p)))) end
      end
    end
    sort!(para,by=x->sum(x.*x)) # increasing "complexity"
# para holds NRPARA random lists of mᵢ,ⱼ defining each a possible A_h
    aA=gcd_partitions(map(p->collectby(i->2*sch[i].mon*p+
           sum(x->x.mon*p*degree(x.pol),sch[i].vcyc),1:length(sch)), para)...)
# aA holds the (a+A)-blocks common to all NRPARA possible A_h
    para=para[1] # choose now the first A_h with "simplest" parameters
    c=map(s->s.coeff*prod(r->iszero(para*r.mon) ? r.pol(r.coeff) : 1,s.vcyc),
          sch)
# c holds the leading coefficients of the Schur elements of A_h
#
# compute now the p-blocks of A_h for all primes p dividing |W|
    res=map(collect(keys(factor(length(W)))))do p
      local bl, vp, j, x, i, cut
      bl=gcd_partitions(blocks(W,p),aA)
# here bl holds the coarsest partition
# - finer than the p-blocks of W
# - finer than the (a+A)-blocks
      if h==[1,-1] return bl end # for 1-parameter algebras bl is h-blocks
      vp=valuation.(c, p)
      i=1
      while i<=length(bl)
# we examine each "pseudo-block" obtained and apply various rules to refine it
        if length(bl[i]) == 1
          if !iszero(vp[bl[i][1]])
            error("Schur elt of v_$p==",vp[bl[i][1]]," alone in block")
          else i+=1
          end
          continue
        end
# a Schur element s is alone in its block iff vₚ[c[s]]=0
        for j in bl[i]
          if vp[j]==0 && length(bl[i])>1
            push!(bl,[j])
            bl[i]=setdiff(bl[i],[j])
          end
        end
# a pseudo-block of size <=3 is a block
        if length(bl[i])<=3
          i+=1
          continue
        end
        if !iszero(h)
# if B is a pseudo-block of size>3 and C ⊂ B is a 0-block
# and for each s in C we have vₚ(c(s))=vₚ(c_0(s)) then C is a block
          for j in filter(x->issubset(x,bl[i]) && x!=bl[i], bl0)
            if vp[j]==map(s->valuation(s.coeff, p), sch[j])
              bl[i]=setdiff(bl[i], j)
              push!(bl, j)
            end
          end
        end
        if length(bl[i])<=3
          i+=1
          continue
        end
# We cut the remaining pseudo-blocks in p-blocks by ultimate test
# ∑_{φ\in bl}φ(T)/s_φ p-integral
        function cut(bl, para)
          local csch, lsch, p0, ct, ch, msch, l, Ah, getH
          InfoChevie("#I p==",p," h",findfirst(==(h),hplanes),":",h," cut",bl,"\n")
          function getH(para)
            o=map(x->x.order,hyperplane_orbits(W))
            o=map(i->sum(o[1:i-1])+(1:o[i]),1:length(o))
            v=map(x->Mvp(:x).^para[x].*E.(length(x),0:length(x)-1),o)
            hecke(W,Tuple(v)) # algebra A_h
          end
# replace para by smallest multiple such that schur elements rational
          l=toM(vcat(map(x->map(y->y.mon,x.vcyc),sch[bl])...))*para
          para*=lcm(denominator.(l))
          Ah=getH(Int.(para))
          csch=map(s->s.coeff*CycPol(Mvp(:x)^(para*s.mon))*
                   prod(x->subs(x.pol,Pol()^(para*x.mon)),s.vcyc),sch[bl])
# csch holds the Schur elements for characters of Ah in bl
          lsch=lcm(csch).//csch
#         InfoChevie(" Schur:", Stime())
          lsch=map(x->x(Mvp(:x)),lsch)
#         InfoChevie(" Value:", Stime(), " ")
          ct=permutedims(CharTable(Ah).irr[bl,:])
          ct=ct[sortperm(classinfo(W)[:classtext],by=x->-length(x)),:]
          nct=count(r->!any(u->u isa Unknown,r),eachrow(ct))
          if nct!=nconjugacy_classes(W)
           xprintln("\n!! Unreliable computation: CharTable($Ah) partially unknown")
          end
          p0=map(x->[x],1:length(bl))
          for ch in eachrow(ct)
            InfoChevie(".")
            msch = lsch.*ch
#           @show Pol.(msch)
#           l=filter(x->valuation(sum(msch[x]),p)>=0,
            l=filter(x->ispInt(sum(msch[x]),p),
              filter(!isempty,map(x->sort(vcat(x...)),combinations(p0))))
            if !(1:length(bl) in l) error("theory ",l,1:length(bl)) end
            l=filter(x->count(y->issubset(y,x),l)==1,l)
#           @show l, p0
            p0=lcm_partitions(l,p0)
            if length(p0)==1
#             InfoChevie(Stime(), " ok\n")
              return [bl]
            end
          end
# here bl has been non-trivially cut
#         InfoChevie(Stime(),"\n  ->",FormatGAP(map(x->bl[x],p0)),"\n")
          return map(x->bl[x],p0)
        end
        j=cut(bl[i], para)
        bl=vcat(bl[1:i-1],j,bl[i+1:end])
        i+=length(j)
      end
      bl=filter(!isempty,bl)
      sort!(bl)
      bl
    end
# The h-blocks is the finest partition coarser than all p-blocks
    res=lcm_partitions(res...)

    #UT: added
    if names
        res = [charnames(W;namesargs...)[i] for i in res]
    end

    if iszero(h) bl0=res end
    return [h,res]
  end
end

"""
    function rouquier_blocks(H::HeckeAlgebra)

The  Rouquier blocks  of a  1-cyclotomic algebra  H is the finest partition
coarser than h-blocks for all hyperplanes h annihilated by H's parameters.
"""
function rouquier_blocks(H::HeckeAlgebra ; names=false, namesargs...)
  W=H.W
  d=rouquier_blocks(W)
  p=vcat(H.para[sort(unique(simple_reps(W)))]...)
  d=filter(x->!isnothing(scalar(prod(p.^(x[1]*lcm(denominator.(x[1])))))),d)
  res = lcm_partitions(map(x->x[2],d)...)

  #UT: added names
  if names
      res = [ [charnames(W;namesargs...)[i] for i in bl] for bl in res ]
  end

  return res
end

end # of module
