#
# grad_tools.jl
#
# Provide tools to load the instrumental parameters and the 
# calibrated data, and for the calculus of the RHAPSODIE data fidelity term. 
#
# ------------------------------------------------
#
# This file is part of Rhapsodie
#
#
# Copyright (c) 2017-2021 Laurence Denneulin (see LICENCE.md)
#

#------------------------------------------------

const ImageInterpolator{T<:AbstractFloat, K<:Kernel{T}} = TwoDimensionalTransformInterpolator{T,K,K}
const MyKer = LinearInterpolators.CatmullRomSpline(Float64, LinearInterpolators.Flat)
#const MyKer = LinearInterpolators.RectangularSpline(Float64, LinearInterpolators.Flat)
#const MyKer = LinearInterpolators.LinearSpline(Float64, LinearInterpolators.Flat)

struct parameters_table{T<: AbstractFloat}
     cols::NTuple{3,Int64} #Imput size (reconstruction)
     rows::NTuple{2,Int64} #Output size (data)
     dataset_length::Int64
     Nframe::Int64
     Nrot::Int64
     Nangle::Int64
     v::Vector{NTuple{2, NTuple{3, T}}}
     indices::Array{Int64,2}
     center::Array{T,1}
     psf_center::NTuple{2,Array{T,1}}
     epsilon::Vector{NTuple{2,Array{T,1}}} 
     derotang::Vector{T}
end


const Trans_Table = Vector{NTuple{2,AffineTransform2D}}(); #Contains all the affine transform used 
const Parameters = parameters_table[];
get_par()::parameters_table = Parameters[1];
const dataset = data_table[];

const PSF_save = Vector{Array{Float64,2}}();
get_PSF()=PSF_save[1];
const EPSILON_save = Array{Float64,1}(undef,1);
function set_epsilon(epsilon)
    EPSILON_save[1]=epsilon
end
get_epsilon()::Float64=EPSILON_save[1];
const PRECOND_SAVE= Vector{Any}();
U()=PRECOND_SAVE[1];

const MASK_save = Vector{Array{Float64,2}}();
get_MASK()=MASK_save[1];


function Indices(S::Int64,Nrot::Int64,Nframe::Int64)
#S = le nombre total de frame gauche ou droite
#Nrot : Nombre de positions du derotateur
#Nframe : Nombre de frame par positions de la lame demi-onde
	
	Nangle=Int32.(S/(Nframe*4)) #Nombre de rotations de la lame demi-onde
	INDICES=zeros(4,Nframe*Nangle)
	for i=1:4
		ind=repeat(range(0,stop=4*Nframe*(Nangle-1),length=Nangle), inner=Nframe)+(Nframe*i .-mod.(range(1,stop=Nframe*Nangle,length=Nframe*Nangle),Nframe))
		INDICES[i,:]=ind
	end
	INDICES=round.(Int64, INDICES);
	return INDICES
end


function Set_Vi(Indices::AbstractArray{Int64,2}; alpha=[0, pi/4, pi/8, 3*pi/8], psi=[0,pi/2])
	J(a)=[cos(2*a) sin(2*a); sin(2*a) -cos(2*a)]
	P(a)=[cos(a), -sin(a)]
    v=Vector{NTuple{2, NTuple{3, Float64}}}(undef, length(Indices));
	for i=1:4
		v1=J(alpha[i])*P(psi[1])
		vn1=norm(v1)^2;
		v2=J(alpha[i])*P(psi[2])
		vn2=norm(v2)^2;
		for k=1:length(Indices[i,:])		
		v[Indices[i,k]] =((round(vn1/2, digits=4), round((abs(v1[1])^2-abs(v1[2])^2)/2, digits=4), round(real(v1[1]*v1[2]), digits=4)),
		                  (round(vn2/2, digits=4), round((abs(v2[1])^2-abs(v2[2])^2)/2, digits=4), round(real(v2[1]*v2[2]), digits=4)));
        end
	end
	return v
end

function Set_Vi(Nframe::Int, dataset_length::Int, mueller_instru::AbstractArray{Float64,2})
    @assert dataset_length == Nframe * size(mueller_instru)[1]
    
    v=Vector{NTuple{2, NTuple{3, Float64}}}(undef, dataset_length);
    for k=1:size(mueller_instru)[1]
        for l=1:Nframe
            v[Nframe*(k-1) + l] =((mueller_instru[k,1], mueller_instru[k,2], mueller_instru[k,3]),
                                  (mueller_instru[k,4], mueller_instru[k,5], mueller_instru[k,6]));
        end
    end
	return v
end


function reset_instrument(V::Vector{NTuple{2, NTuple{3, T}}})  where {T <: AbstractFloat}
    push!(Parameters, parameters_table(get_par().cols, 
                                       get_par().rows, 
                                       get_par().dataset_length, 
                                       get_par().Nframe, 
                                       get_par().Nrot, 
                                       get_par().Nangle, 
                                       V, 
                                       get_par().indices, 
                                       get_par().center, 
                                       get_par().psf_center,
                                       get_par().epsilon, 
                                       get_par().derotang)); 

    popfirst!(Parameters);
end



function TransRotate(A::AffineTransform2D{Float64}, EPSILON, ANGLE, CENTER, NEWCENTER)
	return translate(rotate(translate(CENTER[1]-EPSILON[1],CENTER[2]-EPSILON[2], A), ANGLE),-NEWCENTER[1], -NEWCENTER[2])
end

function bbox_size(inp_dims::NTuple{2,Integer},
                  A::AffineTransform2D{Float64})
                    
	xmin=typemax(Float64);
	ymin=typemax(Float64);
	xmax=typemin(Float64);
	ymax=typemin(Float64);
	
	xmin_int=typemin(Int64);
	ymin_int=typemin(Int64);
	xmax_int=typemin(Int64);
	ymax_int=typemin(Int64);
	
	width=typemin(Float64);
	height=typemin(Float64);
	
	Ind=[repeat(1:inp_dims[1], inner=inp_dims[2]) repeat(1:inp_dims[2], outer=inp_dims[1])]
	for it=1:(inp_dims[1]*inp_dims[2])
		(x,y)=A(Ind[it,1], Ind[it,2]);
		xmin=min(x, xmin);
		ymin=min(y, ymin);
		xmax=max(x, xmax);
		ymax=max(y, ymax);
	end
	xmin_int=floor(Int64,xmin);
	ymin_int=floor(Int64,ymin);
	xmax_int=ceil(Int64,xmax);
	ymax_int=ceil(Int64,ymax);

	width=xmax_int-xmin_int+1+50;
	height=ymax_int-ymin_int+1+50;
    out_dims=(width, height);

	return(out_dims, (xmin_int, xmax_int, ymin_int, ymax_int))
end

function set_fft_op(PSF::AbstractArray{T,2}, PSFCenter::AbstractArray{T,1}) where {T <: AbstractFloat}
 	MapSize=get_par().cols[1:2];
	MapCenter=floor.(MapSize./2).+1
	MAP=zeros(MapSize);
	ker = LinearInterpolators.CatmullRomSpline(Float64, LinearInterpolators.Flat)

	Id = AffineTransform2D{Float64}()
	Recentering=translate(-(MapCenter[1]-PSFCenter[1]), -(MapCenter[2]-PSFCenter[2]), Id)

	LazyAlgebra.apply!(MAP, ker, Recentering, PSF);
	MAP./=sum(MAP);
	push!(PSF_save,MAP);
	F=FFTOperator(MAP)
	FFT=F\Diag(F*ifftshift(MAP)) .*F;
	
	return FFT
end

      

  
#--------------------------------------------UTILS------------------------------------------------
      
function SetCropOperator()
    DSIZE=get_par().rows[1]
    MASK=ones(get_par().cols[1:2]);
    for k=1:length(Trans_Table)
        X11=Trans_Table[k][1](1,1);
        X12=Trans_Table[k][1](1,DSIZE);
        X13=Trans_Table[k][1](DSIZE, DSIZE);
        X14=Trans_Table[k][1](DSIZE, 1);
        
        X21=Trans_Table[k][2](1,1);
        X22=Trans_Table[k][2](1,DSIZE);
        X23=Trans_Table[k][2](DSIZE, DSIZE);
        X24=Trans_Table[k][2](DSIZE, 1);
        Mask1=zeros(get_par().cols[1:2]);
        Mask2=zeros(get_par().cols[1:2]);

        for i=1:get_par().cols[1]
	        for j=1:get_par().cols[2]	
	                H11 = ((i-X11[1])*(X12[1]-X11[1]) +(j-X11[2])*(X12[2]-X11[2]) >0)
	                H12 = ((i-X12[1])*(X13[1]-X12[1]) +(j-X12[2])*(X13[2]-X12[2]) >0)
	                H13 = ((i-X13[1])*(X14[1]-X13[1]) +(j-X13[2])*(X14[2]-X13[2]) >0)
	                H14 = ((i-X14[1])*(X11[1]-X14[1]) +(j-X14[2])*(X11[2]-X14[2]) >0)	
	                H21 = ((i-X21[1])*(X22[1]-X21[1]) +(j-X21[2])*(X22[2]-X21[2]) >0)
	                H22 = ((i-X22[1])*(X23[1]-X22[1]) +(j-X22[2])*(X23[2]-X22[2]) >0)
	                H23 = ((i-X23[1])*(X24[1]-X23[1]) +(j-X23[2])*(X24[2]-X23[2]) >0)
	                H24 = ((i-X24[1])*(X21[1]-X24[1]) +(j-X24[2])*(X21[2]-X24[2]) >0)	
		        if H11 && H12 && H13 && H14
		            Mask1[i,j,:] .=1;
		            #MASK[i,j,:] .=1;
		        end
		        if H21 && H22 && H23 && H24
		            Mask2[i,j,:] .=1;
		            #MASK[i,j,:] .=1;
		        end
	        end
        end
    MASK .*= Mask1.*Mask2
    end
    #write(FitsFile, "MASK.fits", MASK, overwrite=true)
    push!(MASK_save, MASK)
end     

        
function crop(X::M)  where {T<:AbstractFloat, M<:AbstractArray{T,2}}
    #@assert size(X) .==   get_par().cols
    xymin_1=findfirst(get_MASK() .!= 0.)
    xymin_2=findfirst(get_MASK()' .!= 0.)
    xymax_1=findlast(get_MASK() .!= 0.)
    xymax_2=findlast(get_MASK()' .!= 0.)
    
    xmin=min(xymin_1[1], xymin_2[2]);
    ymin=min(xymin_1[2], xymin_2[1]);
    xmax=max(xymax_1[1], xymax_2[2]);
    ymax=max(xymax_1[2], xymax_2[1]);
    
    Y=X[xmin:xmax, ymin:ymax];
    Y[.!isfinite.(Y)].=0    
    return Y     
end        

function crop!(X::M)  where {T<:AbstractFloat, M<:AbstractArray{T,2}}
    #@assert size(X) .==   get_par().cols 
    X[.!isfinite.(X)].=0   
    X .*= get_MASK()        
end        
        
function crop(X::M)  where {T<:AbstractFloat, M<:AbstractArray{T,3}}
    #@assert size(X) .==   get_par().cols
    Y=copy(X);
    return crop!(Y)    
end        

function crop!(X::M)  where {T<:AbstractFloat, M<:AbstractArray{T,3}}
    #@assert size(X) .==   get_par().cols
    for k in size(X)[3] 
        crop!(view(X,:,:,k))
    end     
end        
   
                
function crop(X::PolarimetricMap{T}) where {T<:AbstractFloat}  
    #@assert size(X) .==   get_par().cols
    return PolarimetricMap(X.parameter_type,
                           crop(view(X.I,:,:)),
                           crop(view(X.Q,:,:)),
                           crop(view(X.U,:,:)),
                           crop(view(X.Iu,:,:)),
                           crop(view(X.Ip,:,:)),        
                           crop(view(X.θ,:,:)))   
end        

function crop!(X::PolarimetricMap{T})  where {T<:AbstractFloat}
        crop!(view(X.I,:,:));
        crop!(view(X.Q,:,:));
        crop!(view(X.U,:,:));
        crop!(view(X.Iu,:,:));
        crop!(view(X.Ip,:,:));        
        crop!(view(X.θ,:,:));
end        
        
function pad(X::M)  where {T<:AbstractFloat, M<:AbstractArray{T,2}}
#    X_size = size(X);
#    center_diff = X_size./2 .- get_par().cols[1:2]./2;
#    ker = LinearInterpolators.CatmullRomSpline(Float64, LinearInterpolators.Flat)    
#    Id = AffineTransform2D{Float64}()
#    center_change = translate(center_diff[1],center_diff[2], Id)   
#    PAD=TwoDimensionalTransformInterpolator(get_par().cols[1:2], X_size, ker, ker, center_change)
    PAD = ZeroPaddingOperator(get_par().cols[1:2], size(X)) 
    return PAD*X;
end         
        
function pad(X::PolarimetricMap{T}) where {T<:AbstractFloat}  
    #@assert size(X) .==   get_par().cols
    return PolarimetricMap(X.parameter_type,
                           pad(view(X.I,:,:)),
                           pad(view(X.Q,:,:)),
                           pad(view(X.U,:,:)),
                           pad(view(X.Iu,:,:)),
                           pad(view(X.Ip,:,:)),        
                           pad(view(X.θ,:,:)))   
end    

