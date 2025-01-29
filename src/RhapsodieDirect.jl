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
        load_parameters,
        get_par,
        Load_Data,
        fg!,
        SetCropOperator,
        crop,
        crop!,
        pad,
        set_fft_op,
        TwoDimensionalTransformInterpolator,
        FieldTransformOperator,
        data_simulator,
        generate_model,
        data_generator,
        generate_parameters,
        Double_Difference,
        Double_Ratio,
        Linear_Method,
        NonLinear_Method,
        apply_rhapsodie,
        apply_edge_preserving_smoothing!

    import Base: +, -, *, /, ==, getindex, setindex!, read, write, convert

    using SpecialFunctions
    using TwoDimensional
    using FFTW
    using LinearInterpolators
    using Statistics
    using LinearAlgebra
    using LazyAlgebra
    import LazyAlgebra: Mapping, vcreate, vcopy, apply!
    using StaticArrays
    using FITSIO
    using EasyFITS
    using DelimitedFiles
    #using TiPi
    using Random

    include("Polarimetric_Parameters.jl")
    include("grad_tools.jl")
    include("separable_methods.jl")
    include("datasimul_tools.jl")
end

