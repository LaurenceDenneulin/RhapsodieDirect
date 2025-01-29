#
# Rhapsodie.jl
#
# Package for the Reconstruction of High-contrAst Polarized
# SOurces and Deconvolution for cIrcumstellar Environments (Rhapsodie)
#
#----------------------------------------------------------
#
# 
# Copiyright (c) 2017-2021 Laurence Denneulin (see LICENCE.md)
#

module RhapsodieDirect

    export
        PolarimetricPixel,
        PolarimetricMap,
        write,
        read,
        convert,
        load_field_transforms,
        data_simulator,
        direct_model!,
        data_generator,
        generate_parameters,
        get_indices_table,
        set_default_polarisation_coefficients,
        field_transform,
        bbox_size,
        set_fft_operator,
        set_crop_operator,
        crop,
        pad,
        check_MSE,
        ObjectParameters,
        DatasetParameters,
        FieldTransformParameters,
        FieldTransformOperator

    import Base: +, -, *, /, ==, getindex, setindex!, read, write, convert

    #using SpecialFunctions
    using TwoDimensional
    #using FFTW
    using InterpolationKernels
    using LinearInterpolators
    #using Statistics
    #using LinearAlgebra
    using LazyAlgebra
    import LazyAlgebra: Mapping, vcreate, vcopy, apply!
    #using StaticArrays
    using EasyFITS
    #using DelimitedFiles
    #using Random
    
    include("types.jl")
    include("Polarimetric_Parameters.jl")
    include("utils.jl")
    include("loaders.jl")
    include("separable_methods.jl")
    include("datasimul_tools.jl")
end

