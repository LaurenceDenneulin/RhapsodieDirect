# FieldTransformOperator mappings

function vcreate(::Type{LazyAlgebra.Direct}, A::FieldTransformOperator{T},
                 x::AbstractArray{T,3}, scratch::Bool = false) where {T <: AbstractFloat}
    @assert !Base.has_offset_axes(x)
    @assert size(x) == A.cols
    Array{T,2}(undef, A.rows)
end

function vcreate(::Type{LazyAlgebra.Adjoint}, A::FieldTransformOperator{T},
                 x::AbstractArray{T,2}, scratch::Bool = false) where {T <: AbstractFloat}
    @assert !Base.has_offset_axes(x)
    @assert size(x) == A.rows
    Array{T,3}(undef, A.cols)
end


#FIXME : refactor apply! to make α, β usefull (technically dst = α.R*src + β*dst)

function apply!(α::Real,
                ::Type{LazyAlgebra.Direct},
                R::FieldTransformOperator{T},
                src::AbstractArray{T,3},
                scratch::Bool,
                β::Real,
                dst::AbstractArray{T,2}) where {T<:AbstractFloat}
    @assert β==0 && α==1
    @assert size(src) == R.cols
    @assert size(dst) == R.rows
    n = R.rows[2]
    @assert iseven(n)
    fill!(dst,zero(T));
    # Allocating memory FIXME: find a way to calculate fully in place 
    z = zeros(R.cols[1:2]);
    
    #Compute left direct model
    @simd for i=1:length(R.v_l)
        vupdate!(z, R.v_l[i],view(src,:,:,i)) 
    end
    apply!(view(dst,:, 1:(n÷2)),R.H_l,z);
    
    # Reset the array values to 0. (faster than allocating two different arrays)
    vfill!(z,0.)
    
    # Compute right direct model
    @simd for i=1:length(R.v_r)
        vupdate!(z, R.v_r[i],view(src,:,:,i)) 
    end
    apply!(view(dst,:, (n÷2)+1:n),R.H_r,z);
    return dst
end

function apply!(α::Real,
                ::Type{LazyAlgebra.Adjoint},
                R::FieldTransformOperator{T},
                src::AbstractArray{T,2},
                scratch::Bool,
                β::Real,
                dst::AbstractArray{T,3}) where {T<:AbstractFloat}
    @assert β==0 && α==1
    @assert size(src) == R.rows
    @assert size(dst) == R.cols
    n = R.rows[2]
    @assert iseven(n)
    fill!(dst,zero(T));
  
    y = zeros(R.cols[1:2])
    vmul!(y, R.H_l', view(src, :, 1:(n÷2)))
    @simd for i=1:length(R.v_l)
         vupdate!(view(dst,:,:,i), R.v_l[i], y)
    end
    vmul!(y, R.H_r', view(src, :, (n÷2)+1:n))
    @simd for i=1:length(R.v_r)
         vupdate!(view(dst,:,:,i), R.v_r[i], y)
    end

    return dst;
end

# DirectModel mapping

function vcreate(::Type{LazyAlgebra.Direct}, A::DirectModel{T},
                 x::PolarimetricMap{T}, scratch::Bool = false) where {T <: AbstractFloat}
    @assert !Base.has_offset_axes(x)
    @assert size(x) == A.cols
    Array{T,3}(undef, A.rows)
end

function vcreate(::Type{LazyAlgebra.Adjoint}, A::DirectModel{T},
                 x::AbstractArray{T,3}, scratch::Bool = false) where {T <: AbstractFloat}
    @assert !Base.has_offset_axes(x)
    @assert size(x) == A.rows
    PolarimetricMap{T}(A.parameter_type,Array{T,2}(undef, A.cols),
                                        Array{T,2}(undef, A.cols),
                                        Array{T,2}(undef, A.cols),
                                        Array{T,2}(undef, A.cols),
                                        Array{T,2}(undef, A.cols),
                                        Array{T,2}(undef, A.cols))
end


#FIXME : refactor apply! to make α, β usefull (technically dst = α.R*src + β*dst)

function apply!(α::Real,
                ::Type{LazyAlgebra.Direct},
                R::DirectModel{T},
                src::PolarimetricMap{T},
                scratch::Bool,
                β::Real,
                dst::AbstractArray{T,3}) where {T<:AbstractFloat}
    @assert β==0 && α==1
    @assert size(src) == R.cols
    @assert size(dst) == R.rows
    x = zeros(R.cols[1],R.cols[2], length(src))
    @inbounds for (i,map) in enumerate(get_stokes(src))
        if i>1
            setindex!(x,R.A*map,:,:,i)
        else
            setindex!(x,map,:,:,i)
        end
    end
    
    @inbounds for k=1:length(R.TR)	 
        apply!(view(dst,:,:,k),R.TR[k],x)   
	end
    return dst
end

function apply!(α::Real,
                ::Type{LazyAlgebra.Adjoint},
                R::DirectModel{T},
                src::AbstractArray{T,3},
                scratch::Bool,
                β::Real,
                dst::PolarimetricMap{T}) where {T<:AbstractFloat}
    @assert β==0 && α==1
    @assert size(src) == R.rows
    @assert size(dst) == R.cols
    x = zeros(R.cols[1],R.cols[2], length(dst))
    y = zeros(R.cols[1],R.cols[2], length(dst))
    @inbounds for k=1:length(R.TR)	 
        vmul!(y, R.TR[k]', view(src,:,:,k))
        vupdate!(x,1.,y)   
	end
    @inbounds for (i,map) in enumerate(get_stokes(dst))
        if i>1
        apply!(map,R.A', x[:,:,i])
        else
        vcopy!(map,x[:,:,i])
        end
    end
    rebuild("stokes",dst)
    return dst;
end

