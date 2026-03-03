function load_field_transforms(object::ObjectParameters,
                           data::DatasetParameters{T},
                           parameters::Vector{FieldTransformParameters}) where {T<:AbstractFloat}
    @assert data.frames_total == length(parameters)
                           
    Id = AffineTransform2D{T}()
    field_transforms=Vector{FieldTransformOperator{T}}()
    
    for k=1:data.frames_total
        T_left=field_transform(Id, 
                               parameters[k].translation_left, 
                               parameters[k].field_angle, 
                               object.center,data.center) 
        T_right=field_transform(Id, 
                               parameters[k].translation_right, 
                               parameters[k].field_angle, 
                               object.center,data.center)   

                                   
        output_size=(data.size[1], data.size[2]÷2)
        input_size= object.size
    	T1=TwoDimensionalTransformInterpolator(output_size, 
    	                                       input_size, 
    	                                       parameters[k].ker, 
    	                                       parameters[k].ker, 
    	                                       inv(T_left))
    	T2=TwoDimensionalTransformInterpolator(output_size, 
    	                                       input_size, 
    	                                       parameters[k].ker, 
    	                                       parameters[k].ker, 
    	                                       inv(T_right))
    	
	    push!(field_transforms, FieldTransformOperator((object.size[1],object.size[2],3), 
	                                                     data.size, 
	                                                     parameters[k].polarization_left,
	                                                     parameters[k].polarization_right,
	                                                     T1,
	                                                     T2));    

    end
    #TODO: Add bounding_box and mask calculus (cf. bbox_size and SetCropOperator) and export it somehow
    return field_transforms
end        

