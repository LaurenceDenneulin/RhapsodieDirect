function load_field_transforms(object::ObjectParameters,
                           data::DatasetParameters,
                           parameters::Vector{FieldTransformParameters})
    @assert data.frames_total == length(parameters)
                           
    Id = AffineTransform2D{Float64}()
    field_transforms=Vector{FieldTransformOperator{Float64}}()
    
    for k=1:data.frames_total
        T_left=field_transform(Id, 
                               parameters[k].translation_left, 
                               parameters[k].field_angle, 
                               data.center, object.center)
        T_right=field_transform(Id, 
                               parameters[k].translation_right, 
                               parameters[k].field_angle, 
                               data.center, object.center)   

                           
        output_size=(data.size[1], data.size[2]÷2)
        input_size= object.size
    	T1=TwoDimensionalTransformInterpolator(output_size, 
    	                                       input_size, 
    	                                       parameters[k].ker, 
    	                                       parameters[k].ker, 
    	                                       T_left)
    	T2=TwoDimensionalTransformInterpolator(output_size, 
    	                                       input_size, 
    	                                       parameters[k].ker, 
    	                                       parameters[k].ker, 
    	                                       T_right)
    	
	    push!(field_transforms, FieldTransformOperator((object.size[1],object.size[2],3), 
	                                                     data.size, 
	                                                     parameters[k].polarization_left,
	                                                     parameters[k].polarization_right,
	                                                     T1,
	                                                     T2));    

    end
    return field_transforms
end        

