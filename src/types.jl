"""

""" ObjectParameters
struct ObjectParameters{T<:AbstractFloat, U<:Int}
    size::NTuple{2,U}
    center::NTuple{2,T}
end

"""

""" DatasetParameters
struct DatasetParameters{T<:AbstractFloat, U<:Int}
    size::NTuple{2,U}
    frames_total::U
    frames_per_hwp_pos::U
    hwp_cycles::U
    center::NTuple{2,T}
end


"""

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

""" FieldTransformOperator
struct FieldTransformOperator{T<:AbstractFloat, L<:Mapping, R<:Mapping} <: LinearMapping
    cols::NTuple{3,Int} 
    rows::NTuple{2,Int} 
    v_l::NTuple{3,T}
    v_r::NTuple{3,T}
    H_l::L              
    H_r::R      
end

function FieldTransformOperator(cols::NTuple{3,Int64},
                                rows::NTuple{2,Int64},
                                v_l::NTuple{3,T},
                                v_r::NTuple{3,T},
                                T_l::TwoDimensionalTransformInterpolator{T},
                                T_r::TwoDimensionalTransformInterpolator{T}) where {T <: AbstractFloat}
    FieldTransformOperator(cols, rows, v_l, v_r, T_l, T_r)
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

    z_l = zeros(R.cols[1:2]);
    z_r = zeros(R.cols[1:2]);
    for i=1:length(R.v_l)
        # FIXME: use Generalized Matrix-Vector Multiplication of LazyAlgebra here.
        x_i = view(src,:,:,i)
	z_l .+= R.v_l[i]*x_i;
	z_r .+= R.v_r[i]*x_i;
    end
    dst[:, 1:(n÷2)  ] = R.H_l*z_l;
    dst[:, (n÷2)+1:n] = R.H_r*z_r;
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
    y_l = R.H_l'*view(src, :, 1:(n÷2))
    y_r = R.H_r'*view(src, :, (n÷2)+1:n)
    for i=1:length(R.v_l)
        #dst[:,:,i] = R.v_l[i]*y_l + R.v_r[i]*y_r;
        # FIXME: The following should be equivalent and much faster:
         vcombine!(view(dst,:,:,i), R.v_l[i], y_l, R.v_r[i], y_r)
    end
    return dst;
end

