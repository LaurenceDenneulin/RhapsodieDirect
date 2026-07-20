# FieldTransformOperator mappings
function unsafe_vmul!(α::Number,
                R::FieldTransformOperator{T},
                src::AbstractVector{AbstractMatrix{T}},
                β::Number,
                dst::AbstractArray{T,2}) where {T<:AbstractFloat}
    n = R.rows[2]
    @assert iseven(n)
    # Compute left direct model
    z = zeros(T,R.cols[1:2])
    for i=1:length(R.v_l)
        vupdate!(z, α*R.v_l[i],src[i])  
        #FIXME: update H to new API (TwoDimentionalInterpolator)
    end
    vmul!(ONE, R.H_l, z,β, view(dst,:, 1:(n÷2)))

    # Compute right direct model
    vzeros!(z)   
    for i=1:length(R.v_r)
        vupdate!(z, α*R.v_r[i],src[i]) 
    end
    vmul!(ONE, R.H_r, z,β, view(dst,:, (n÷2)+1:n))
    return nothing
end

function unsafe_vmul!(α::Real,
                R::LazyAlgebra.Adjoint{<:FieldTransformOperator{T}},
                src::AbstractArray{T,2},
                β::Real,
                dst::AbstractVector{AbstractMatrix{T}}) where {T<:AbstractFloat}
    pR=parent(R)
    n = pR.rows[2]
    @assert iseven(n)
    y = Array{T}(undef,pR.cols[1:2])
    # Compute left adjoint model     
    vmul!(ONE, pR.H_l', view(src, :, 1:(n÷2)), ZERO, y)
    for i=1:length(pR.v_l)
         vcombine!(α*pR.v_l[i],y, β, dst[i])
    end
    
    # Compute right adjoint model
    vmul!(ONE, pR.H_r', view(src, :, (n÷2)+1:n), ZERO, y)
    for i=1:length(pR.v_r)
         vcombine!(α*pR.v_r[i],y, ONE, dst[i])
    end
    return nothing
end

# LinearDirectModel mapping

#TODO : Modify for new polarimetricmap structure
function unsafe_vmul!(α::Real,
                R::LinearDirectModel{T},
                src::PolarimetricMap{T},
                β::Real,
                dst::AbstractArray{T,3}) where {T<:AbstractFloat}
    x = zeros(T,R.cols[1],R.cols[2], length(src))
    @inbounds for (i,map) in enumerate(get_stokes(src))
        setindex!(x,R.A*map,:,:,i)
    end
    
    @inbounds for k=1:length(R.TR)	 
        vmul!(α,R.TR[k],x,β,view(dst,:,:,k))   
	end
    return nothing
end

function unsafe_vmul!(α::Real,
                R::LazyAlgebra.Adjoint{LinearDirectModel{T}},
                src::AbstractArray{T,3},
                β::Real,
                dst::PolarimetricMap{T}) where {T<:AbstractFloat}
    pR=parent(R)
    x = convert(dst)#zeros(T,R.cols[1],R.cols[2], length(dst))
    vmul!(ONE, pR.TR[1]', view(src,:,:,k),β, x)
    @inbounds for k=2:length(pR.TR)	 
        vmul!(ONE, pR.TR[k]', view(src,:,:,k),ONE,x)
    end
    @inbounds for (i,map) in enumerate(get_stokes(dst))
        vmul!(α,R.A', x[:,:,i], ZERO, map)
    end
    rebuild("stokes",dst)
    return dst
end

