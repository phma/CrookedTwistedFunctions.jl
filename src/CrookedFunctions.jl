module CrookedFunctions
# Rust code by Hitokiri Battosai aka Ferecides de Siros
# Rewritten in Julia by Pierre Abbat
# Public domain.
using OffsetArrays
export derivatives

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

end # module CrookedFunctions
