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
                               
#------------------------------------------------
# Constructors 
 
function PolarimetricMap(parameter_type::AbstractString, 
                         x1::A, 
                         x2::A, 
                         x3::A) where {T<:AbstractFloat, A<:AbstractArray{T,2}}                   
    n1, n2 = size(x1)
    @assert ((n1,n2) == size(x2)) && ((n1,n2) == size(x3)) 
        if parameter_type == "stokes"
        I = x1            # Stokes parameter I (total light intensity)
        Q = x2            # Stokes parameter Q
        U = x3            # Stokes parameter U
        Iu = Array{T}(undef, n1, n2)    # intensity of unpolarized light
        Ip = Array{T}(undef, n1, n2)    # intensity of linearly polarized light
        θ = Array{T}(undef, n1, n2) # angle of linearly polarized light
        @inbounds for i2 in 1:n2
            @simd for i1 in 1:n1
                Ip[i1,i2] = sqrt(Q[i1,i2]^2 + U[i1,i2]^2);
                Iu[i1,i2] = I[i1, i2]-Ip[i1,i2]; 
                θ[i1,i2] = atan(U[i1,i2], Q[i1,i2])/2;
            end
        end
    elseif parameter_type == "intensities"
        I = Array{T}(undef, n1, n2)    # Stokes parameter I (total light intensity)
        Q = Array{T}(undef, n1, n2)    # Stokes parameter Q
        U = Array{T}(undef, n1, n2)    # Stokes parameter U
        Iu = x1          # intensity of unpolarized light
        Ip = x2          # intensity of linearly polarized light
        θ = x3       # angle of linearly polarized light
         @inbounds for i2 in 1:n2
            @simd for i1 in 1:n1
                I[i1,i2] = Iu[i1,i2] + Ip[i1,i2];
                Q[i1,i2] = Ip[i1,i2] .*cos.(2*θ[i1,i2]);
                U[i1,i2] = Ip[i1,i2] .*sin.(2*θ[i1,i2]);
            end
        end
    elseif parameter_type == "mixed" 
        I = Array{T}(undef, n1, n2)     # Stokes parameter I (total light intensity
        Q = x2            # Stokes parameter Q
        U = x3            # Stokes parameter U
        Iu = x1          # intensity of unpolarized light     
        Ip = Array{T}(undef, n1, n2)    # intensity of linearly polarized light
        θ = Array{T}(undef, n1, n2) # angle of linearly polarized light
        @inbounds for i2 in 1:n2
            @simd for i1 in 1:n1
                Ip[i1,i2] = sqrt(Q[i1,i2]^2 + U[i1,i2]^2);
                I[i1,i2] = Iu[i1,i2] + Ip[i1,i2];
                θ[i1,i2] = atan(U[i1,i2], Q[i1,i2])/2;
            end
        end

    else
        error("unkown type, only known types : stokes, intensities and mixed")
    end
    PolarimetricMap{T}(parameter_type,I,Q,U,Iu,Ip,θ)
end

function PolarimetricMap(parameter_type::AbstractString, 
                         x::Array{T,3}) where {T<:AbstractFloat}
    n1, n2, n3 = size(x)
    @assert n3 == 3 
    PolarimetricMap(parameter_type, 
                    copy(x[:,:,1]),
                    copy(x[:,:,2]),
                    copy(x[:,:,3]));
end

function PolarimetricMap{T}(parameter_type::AbstractString, n1::Int, n2::Int) where {T<:AbstractFloat}
    return PolarimetricMap(parameter_type,
                           Array{T,2}(undef, n1, n2),
                           Array{T,2}(undef, n1, n2),
                           Array{T,2}(undef, n1, n2),
                           Array{T,2}(undef, n1, n2),
                           Array{T,2}(undef, n1, n2),
                           Array{T,2}(undef, n1, n2))
end

function (P::PolarimetricMap)(X::PolarimetricMap)
    P.I .= copy(X.I)
    P.Q .= copy(X.Q)
    P.U .= copy(X.U)
    P.Iu .= copy(X.Iu)
    P.Ip .= copy(X.Ip)
    P.θ .= copy(X.θ)
end


"""
    
""" rebuild
function rebuild(parameter_type::AbstractString,P::PolarimetricMap)
    if parameter_type == "stokes"
        P.Ip .= sqrt.(P.Q.^2 +P.U.^2)
        P.Iu .= P.I - P.Ip
        P.θ .= atan.(P.U,P.Q)/2          
    elseif parameter_type == "intensities"
        P.I .= P.Iu + P.Ip
        P.Q .= P.Ip .* cos.(2*P.θ)
        P.U .= P.Ip .* sin.(2*P.θ)
    elseif parameter_type == "mixed" 
        P.Ip .= sqrt.(P.Q.^2 +P.U.^2)
        P.I .= P.Iu + P.Ip
        P.θ .= atan.(P.U,P.Q)/2 
    else
        error("unkown type, only known types : stokes, intensities and mixed")
    end
    #(minimum(P.I) < 0) && @warn "Negative intensities in I"  
    #(minimum(P.Iu) < 0) && @warn "Negative intensities in Iu"  
    #(minimum(P.Ip) < 0) && @warn "Negative intensities in Ip"  
    return nothing
end


function (P::PolarimetricMap)(X::AbstractArray{T,3}) where {T<: AbstractFloat}
    @inbounds for (i,map) in enumerate(P)
        map .= copy(X[:,:,i])
    end
    rebuild(P.parameter_type,P)
end

#------------------------------------------------
# Base fonction redefinitions

Base.size(A::PolarimetricMap) = size(A.I)
Base.length(A::PolarimetricMap) = 3
Base.copy(X::PolarimetricMap{T}) where {T<:AbstractFloat} = PolarimetricMap(X.parameter_type, 
                                                         copy(X.I), 
                                                         copy(X.Q), 
                                                         copy(X.U), 
                                                         copy(X.Iu), 
                                                         copy(X.Ip), 
                                                         copy(X.θ))

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


function read(T::Type{ST},parameter_type::AbstractString, filename::AbstractString) where {ST<:AbstractFloat}
    X=readfits(filename);
    return PolarimetricMap{T}(parameter_type, 
                           view(X,:,:,4)', 
                           view(X,:,:,5)', 
                           view(X,:,:,6)', 
                           view(X,:,:,1)', 
                           view(X,:,:,2)', 
                           view(X,:,:,3)')
end

read(parameter_type::AbstractString, 
     filename::AbstractString) = read(Float64,
                                      parameter_type::AbstractString, 
                                      filename::AbstractString)


