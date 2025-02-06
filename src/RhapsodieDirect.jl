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
        convert,
        data_generator,
        data_simulator,
        DatasetParameters,
        direct_model!,
        field_transform,
        FieldTransformOperator,
        FieldTransformParameters,
        generate_parameters,
        get_indices_table,
        load_field_transforms,
        ObjectParameters,
        PolarimetricMap,
        PolarimetricPixel,
        read,
        set_default_polarisation_coefficients,
        set_fft_operator,
        write

    import Base: +, -, *, /, ==, getindex, setindex!, read, write, convert

    using TwoDimensional
    using FFTW
    using InterpolationKernels
    using LinearInterpolators
    using LazyAlgebra
    import LazyAlgebra: Mapping, vcreate, vcopy, apply!
    using EasyFITS
    
    include("types.jl")
    include("polarimetric_parameters.jl")
    include("mappings.jl")
    include("model.jl")
    include("utils.jl")
    include("loaders.jl")
    include("datasimul_tools.jl")
end

