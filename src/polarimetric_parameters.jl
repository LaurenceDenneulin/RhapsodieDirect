#
# Polarimetric_Parameters.jl
#
# Provide the type polarimetric parameter.
#
# ------------------------------------------------
#
# This file is part of Rhapsodie
#
#
# Copyright (c) 2017-2021 Laurence Denneulin (see LICENCE.md)

# List for parameters order depending of the type
const ParameterTypes = ((:I,:Q,:U), #Stokes (linear)
                        (:Iu,:Q,:U), #Mixed (not linear)
                        (:Iu,:Ip,:θ), #Intensities (not linear)
                        (:Iu,:Ip), #Intensities with theta fixed (linear)
                        (:I_star,:I_disk, :Q,:U), #Stokes (linear) with ADI
                        (:Iu_star,:Iu_disk, :Q,:U), #Mixed (not linear) with ADI
                        (:Iu_star,:Iu_disk,:Ip,:θ), #Intensities (not linear) with ADI
                        (:Iu_star,:Iu_disk,:Ip)) #Intensities with theta fixed (linear) with ADI
#------------------------------------------------
"""
PolarimetricMap{S}(x) -> PolarimetricMap

create an object of type PolarimetricMap where:
    - x is a vector of matrix
    - S is a Tuple of Symbol indicating the type of map as in the following list : $(ParameterTypes)

For exemple with a for Stokes parameters S=(I,Q,U) given as a vector of matrices:

using Rhapsodie
X = PolarimetricMap{(:I,:Q,:U)}(S)
X.I #yields the Stokes parameter I
X.Ip #yields the polarized intensity Ip
X.I[1,1] #yields the Stokes parameter I at the CartesianIndex (1,1); 

""" PolarimetricMap   
struct PolarimetricMap{S, T<: AbstractFloat, V<:AbstractVector{<:AbstractMatrix{T}}}
    data::V
    function PolarimetricMap{S}(data::AbstractVector{<:AbstractMatrix{T}}) where {S,T<:AbstractFloat}
        S isa Tuple{Vararg{Symbol}} || error("S must be a tuple of symbols")
        axes(data,1) == 1:length(S) || error("data and field names must have the same indices")
        shape = axes(data[1])
        for d in 2:length(S)
            axes(data[d]) == shape || error("all data fields must have the same shape")
        end
        return new{S,T,typeof(data)}(data)
    end
end

function PolarimetricMap{S}(data::M) where {S,T<:AbstractFloat, M<:AbstractArray{T,3}}
    data_vect=Vector{Matrix{T}}()
    for i=1:axes(data,3)
        push!(data_vect, data[:,:,i])
    end
    PolarimetricMap{S}(data_vect);
end


@generated function Base.propertynames(x::PolarimetricMap{S}) where {S}
    return quote
        $(Tuple(sort([S..., :data])))
    end
end

@inline Base.getproperty(x::PolarimetricMap, key::Symbol) = _getproperty(x, Val(key))

_getproperty(x::PolarimetricMap, ::Val{:data}) = getfield(x, :data)
_getproperty(x::PolarimetricMap, ::Val{key}) where {key} = throw(KeyError(key))

# trait

polarimetricfields(x::PolarimetricMap) = polarimetricfields(typeof(x))
polarimetricfields(::Type{<:PolarimetricMap{S,T,A}}) where {S,T,A} = S

#Base.eltype(x::PolarimetricMap) = polarimetricfields(typeof(x))
Base.eltype(::Type{<:PolarimetricMap{S,T,A}}) where {S,T,A} = T

for S in ParameterTypes
    for i in 1:length(S)
        @eval begin
            _getproperty(x::PolarimetricMap{$S,T}, ::$(Val{S[i]})) where {T} = getfield(x, :data)[$i]
        end
    end
end
   
_getproperty(x::PolarimetricMap{(:I,:Q,:U)}, ::Val{:Ip}) = sqrt.(x.Q.^2 + x.U.^2)   
_getproperty(x::PolarimetricMap{(:I,:Q,:U)}, ::Val{:Iu}) = x.I - x.Ip
_getproperty(x::PolarimetricMap{(:I,:Q,:U)}, ::Val{:θ}) = atan.(x.U, x.Q)   

_getproperty(x::PolarimetricMap{(:Iu,:Ip,:θ)}, ::Val{:I}) = x.Iu + x.Ip
_getproperty(x::PolarimetricMap{(:Iu,:Ip,:θ)}, ::Val{:Q}) = x.Ip.*cos.(2*x.θ)
_getproperty(x::PolarimetricMap{(:Iu,:Ip,:θ)}, ::Val{:U}) = x.Ip.*sin.(2*x.θ)
 
_getproperty(x::PolarimetricMap{(:I_star,:I_disk,:Q,:U)}, ::Val{:Ip}) = sqrt.(x.Q.^2 + x.U.^2)
_getproperty(x::PolarimetricMap{(:I_star,:I_disk,:Q,:U)}, ::Val{:I}) = x.I_star + x.I_disk  
_getproperty(x::PolarimetricMap{(:I_star,:I_disk,:Q,:U)}, ::Val{:Iu_star}) = x.I_star    
_getproperty(x::PolarimetricMap{(:I_star,:I_disk,:Q,:U)}, ::Val{:Iu}) = x.I - x.Ip 
_getproperty(x::PolarimetricMap{(:I_star,:I_disk,:Q,:U)}, ::Val{:Iu_disk}) = x.I_disk - x.Ip
_getproperty(x::PolarimetricMap{(:I_star,:I_disk,:Q,:U)}, ::Val{:θ}) = atan.(x.U, x.Q)   

_getproperty(x::PolarimetricMap{(:Iu_star,:Iu_disk,:Ip,:θ)}, ::Val{:Iu}) = x.Iu_star + x.Iu_disk
_getproperty(x::PolarimetricMap{(:Iu_star,:Iu_disk,:Ip,:θ)}, ::Val{:I}) = x.Iu + x.Ip
_getproperty(x::PolarimetricMap{(:Iu_star,:Iu_disk,:Ip,:θ)}, ::Val{:I_star}) = x.Iu_star
_getproperty(x::PolarimetricMap{(:Iu_star,:Iu_disk,:Ip,:θ)}, ::Val{:I_disk}) = x.Iu_disk + x.Ip
_getproperty(x::PolarimetricMap{(:Iu_star,:Iu_disk,:Ip,:θ)}, ::Val{:Q}) = x.Ip.*cos.(2*x.θ)
_getproperty(x::PolarimetricMap{(:Iu_star,:Iu_disk,:Ip,:θ)}, ::Val{:U}) = x.Ip.*sin.(2*x.θ)
   
                               
#TODO: finish this function if interesting.
#=
function (P::PolarimetricMap{S1})(X::PolarimetricMap{S2})
    @assert S1 == S2 || error("Polarimetric parameters must be of the same type")
    P.data .= copy(X.data)
    
end
=#

#------------------------------------------------
# Base fonction redefinitions

Base.size(A::PolarimetricMap) = size(A[1])
Base.length(A::PolarimetricMap{S}) where {S} = length(S)
Base.copy(X::PolarimetricMap{S}) where {S} = PolarimetricMap{S}(X.data)
#=
function get(x::PolarimetricMap{T}, i::Int64) where {T<:AbstractFloat}
    return eval(:($(x).$(MAPDICT[x.parameter_type][i])))
end

function Base.iterate(x::PolarimetricMap{T}, state=1) where {T<:AbstractFloat}
    n = length(x)
    if state > length(x)
        return nothing
    else
        return (get(x,state),state+1)
    end
end

function get_stokes(x::PolarimetricMap{T}) where {T<:AbstractFloat}
    return [x.I,x.Q,x.U]
end

function get_mixed(x::PolarimetricMap{T}) where {T<:AbstractFloat}
    return [x.Iu,x.Q,x.U]
end

function get_intensities(x::PolarimetricMap{T}) where {T<:AbstractFloat}
    return [x.Iu,x.Ip,x.θ]
end

function Base.fill!(P::PolarimetricMap, r::Real)
    @inbounds for (i,map) in enumerate(P)
        map .= r
    end
    rebuild(P.parameter_type,P)
end

 function +(x::PolarimetricMap, y::Array{T,3}) where {T<:AbstractFloat} 
    @assert size(y)[1:2] == size(x)       
    if x.parameter_type == "stokes"
       I=x.I + view(y,:,:,1);
       Q=x.Q + view(y,:,:,2);
       U=x.U + view(y,:,:,3);
       return PolarimetricMap("stokes", I, Q, U)
    elseif x.parameter_type == "intensities"
       Iu=x.Iu + view(y,:,:,1);
       Ip=x.Ip + view(y,:,:,2);
       θ=x.θ + view(y,:,:,3);
       return PolarimetricMap("intensities", Iu, Ip, θ)
    elseif x.parameter_type == "mixed"
       Iu=x.Iu + view(y,:,:,1);
       Q=x.Q + view(y,:,:,2);
       U=x.U + view(y,:,:,3);
       return PolarimetricMap("mixed", Iu, Q, U)
    else
        error("unknown parameter type")
    end
 end
 
 +(y::Array{T,3}, x::PolarimetricMap) where {T<:AbstractFloat} = x + y
 -(x::PolarimetricMap, y::Array{T,3}) where {T<:AbstractFloat} = x + (-y)
 

function +(x::PolarimetricMap{T1}, y::PolarimetricMap{T2}) where {T1<:AbstractFloat,T2<:AbstractFloat} 
    if x.parameter_type != y.parameter_type
        @warn "x.parameter_type : "*x.parameter_type*" is different of y.parameter_type : "*y.parameter_type*". The result of the sum will be of parameter_type : "*x.parameter_type*"."
    end
    return x + convert(Array{T1,3}, y, x.parameter_type)
end                  


 function -(x::PolarimetricMap{T1}, y::PolarimetricMap{T2}) where {T1<:AbstractFloat,T2<:AbstractFloat} 
    if x.parameter_type != y.parameter_type
        @warn "x.parameter_type : "*x.parameter_type*" is different of y.parameter_type : "*y.parameter_type*". The result of the sum will be of parameter_type : "*x.parameter_type*"."
    end
    return x - convert(Array{T1,3}, y, x.parameter_type)
end
  
 vcopy(x::PolarimetricMap) = PolarimetricMap(x.parameter_type,
                                             x.I, 
                                             x.Q, 
                                             x.U,
                                             x.Iu,
                                             x.Ip, 
                                             x.θ)
 function vcreate(x::PolarimetricMap)
 
    @assert (x.parameter_type == "stokes") | 
            (x.parameter_type == "intensities") | 
            (x.parameter_type == "mixed")
     n1,n2=size(x);
     return PolarimetricMap(x.parameter_type, n1, n2)
 end     
                                    
 function convert(::Type{Array{T,3}}, x::PolarimetricMap{T}, parameter_type::AbstractString) where {T <:AbstractFloat}
     if parameter_type == "stokes"
       return cat(x.I, x.Q, x.U, dims=3)
    elseif parameter_type == "intensities"
       return cat(x.Iu, x.Ip, x.θ, dims=3)
    elseif parameter_type == "mixed"
       return cat(x.Iu, x.Q, x.U, dims=3)
    else
        error("unknown parameter type")
    end
 end
 
 function convert(::Type{Array{T,3}}, x::PolarimetricMap{T}) where {T <:AbstractFloat}
     convert(Array{T,3},x, x.parameter_type)
 end

function vnorm2(x::PolarimetricMap{T}) where {T <:AbstractFloat}
    @assert (x.parameter_type == "stokes") | 
            (x.parameter_type == "intensities") | 
            (x.parameter_type == "mixed")
     n1,n2=size(x);
    if x.parameter_type == "stokes"
        return (vdot(x.I,x.I) + vdot(x.Q,x.Q) + vdot(x.U,x.U))/3
    elseif x.parameter_type == "intensities"
        return (vdot(x.Iu,x.Iu) + vdot(x.Ip,x.Ip) + vdot(x.θ,x.θ))/3
    elseif x.parameter_type == "mixed"
        return (vdot(x.Iu,x.Iu) + vdot(x.Q,x.Q) + vdot(x.U,x.U))/3
    end
end
=#

#------------------------------------------------
# Writting function to save PolarimetricMap in fits file
"""
    write(X,'filename.fits') 
    
where X is a PolarimetricMap, write a fitsfile

"""

function write(X::PolarimetricMap, filename::AbstractString)
    data=cat(X.Iu', X.Ip', X.θ', X.I', X.Q', X.U',dims=3)
 
    writefits!(filename,
              ["MAPORDER" => "Iu, Ip, Theta, I, Q, U"],
              data)
end

"""
    read('parameter_type','filename.fits') -> PolarimetricMap
    
create an object of type PolarimetricMap from a fits file with:
    - Parameters I, Q, U (i.e. parameter_type = 'stokes')
    - Parameters Iu, Ip and θ (i.e. parameter_type = 'intensities')
    - Parameters Iu, Q, U (i.e. parameter_type = 'mixed')
   
"""


function read(T::Type{ST},S, filename::AbstractString) where {ST<:AbstractFloat,S}
    X=readfits(filename);
    return PolarimetricMap{S}(X)
end

read(parameter_type::AbstractString, 
     filename::AbstractString) = read(Float64,
                                      parameter_type::AbstractString, 
                                      filename::AbstractString)


