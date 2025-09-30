module CrookedFunctions
# Rust code by Hitokiri Battosai aka Ferecides de Siros
# Rewritten in Julia by Pierre Abbat
# Public domain.
using OffsetArrays
export derivatives,isPermutation,twist

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

function twist(order::Integer)
  sz=2^order
  ret=OffsetVector(fill(0,sz),-1)
  for i in eachindex(ret)
    cnt=count_ones(i)
    ret[i]=(i<<cnt|i>>(order-cnt))&(sz-1)
  end
  ret
end

end # module CrookedFunctions
