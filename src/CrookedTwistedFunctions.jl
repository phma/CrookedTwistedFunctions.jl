module CrookedTwistedFunctions
# Rust code by Hitokiri Battosai aka Ferecides de Siros
# Rewritten in Julia by Pierre Abbat
# Public domain.
using OffsetArrays
export derivatives,isPermutation,apnScore,twist

"""
    derivatives(sbox::OffsetVector{<:Integer})

Given a vector of 2^n integers starting at index 0, compute the number of
different values the derivative takes for each nonzero difference.
"""
function derivatives(sbox::OffsetVector{<:Integer})
  sz=length(sbox)
  derivativeCounts=Int[]
  for a in 1:sz-1
    der=Set{eltype(sbox)}()
    for x in 0:sz-1
      push!(der,sbox[x⊻a]⊻sbox[x])
    end
    push!(derivativeCounts,length(der))
  end
  derivativeCounts
end

"""
    isPermutation(sbox::OffsetVector{<:Integer})

Return true if the vector is a permutation of its indices.
"""
function isPermutation(sbox::OffsetVector{<:Integer})
  # isperm expects numbers starting at 1
  sorted=sort(sbox)
  for i in eachindex(sorted)
    if sorted[i]!=i
      return false
    end
  end
  return true
end

"""
    apnScore(derivativeCounts::Vector{<:Integer})

Compute the APN (almost perfectly nonlinear) score of an S-box, given the output
of `derivatives`.
"""
function apnScore(derivativeCounts::Vector{<:Integer})
  expected=(length(derivativeCounts)+1)÷2
  perfectMatches=count(x->x==expected,derivativeCounts)
  perfectMatches/length(derivativeCounts)
end

"""
    apnScore(derivativeCounts::Vector{<:Integer})

Compute the APN (almost perfectly nonlinear) score of an S-box, whose domain and
range should be the same size.
"""
function apnScore(sbox::OffsetVector{<:Integer})
  apnScore(derivatives(sbox))
end

"""
    twist(order::Integer)

Compute the canonical twisted function of order `order` as a 0-based vector.
The canonical twisted function is the one with no exclusive-or or bit permutation
steps, just rotating the bit vector by the number of one-bits.
"""
function twist(order::Integer)
  sz=2^order
  ret=OffsetVector(fill(0,sz),-1)
  for i in eachindex(ret)
    cnt=count_ones(i)
    ret[i]=(i<<cnt|i>>(order-cnt))&(sz-1)
  end
  ret
end

end # module CrookedTwistedFunctions
