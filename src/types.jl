const ImageInterpolator{T<:AbstractFloat, K<:Kernel{T}} = TwoDimensionalTransformInterpolator{T,K,K}

"""
    ObjectParameters(size, center)
    
where:

* `size` is the size of one map of the object
* `center` is the center of the object. 

""" ObjectParameters
struct ObjectParameters{T<:AbstractFloat, U<:Int}
    size::NTuple{2,U}
    center::NTuple{2,T}
end

"""
    DatasetParameters(size, frames_total, frames_per_hwp_pos, hwp_cycles, center)
    
* `size` is the total size of one frame on ne camera (with both side)
* `frames_total` is the total number of frames in the dataset
* `frames_per_hwp_pos`is the number of frames for a given half-wave-plate (hwp) position
* `hwp_cycles` is the number of hwp cycles
* `center` is the center of the object in this dataset

""" DatasetParameters
struct DatasetParameters{T<:AbstractFloat, U<:Int}
    size::NTuple{2,U}
    frames_total::U
    frames_per_hwp_pos::U
    hwp_cycles::U
    center::NTuple{2,T}
end


"""
    FieldTransformParameters(ker,field_angle, translation_left, translation_right, polarization_left, polarization_right)

* `ker` is an interpolation Kernel from the package `InterpolationKernels`
* `field_angle` is the field rotation for the given frame
* `translation_left` is the vector of translation composed with the difference of the object_center on the left side of the camera and the dataset center.
* `translation_right`is the vector of translation composed with the difference of the object_center on the right side of the camera and the dataset center.
* `polarization_left` are the polarization coefficient (the three first mueller matrix coefficients) of the left side of the camera
* `polarization_right`are the polarization coefficient (the three first mueller matrix coefficients) of the right side of the camera

""" FieldTransformParameters
struct FieldTransformParameters{T<:AbstractFloat,K<:Kernel}
    ker::K
    field_angle::T
    translation_left::NTuple{2,T}
    translation_right::NTuple{2,T}
    polarization_left::NTuple{3,T}
    polarization_right::NTuple{3,T}
end


"""
    FieldTransformOperator
    
provides the linear combination of the geometrical transform and polarization coefficient, with the blurred Stokes parameters.

""" FieldTransformOperator
struct FieldTransformOperator{T<:AbstractFloat, 
                              ColType<:NTuple{3,Int},
                              RowType<:NTuple{2,Int},
                              P<:NTuple{3,T},
                              L<:Mapping, 
                              R<:Mapping} <: LinearMapping
    cols::ColType
    rows::RowType
    v_l::P
    v_r::P
    H_l::L              
    H_r::R      
end

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

