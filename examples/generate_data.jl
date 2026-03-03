using Revise
using RhapsodieDirect
using DelimitedFiles
using AstroFITS
using InterpolationKernels

if prod(readdir() .!= "test_results")     
    mkdir("test_results")
end
T=Float32

ker=CatmullRomSpline(T, Flat)
tau=T(0.25)
object_params=ObjectParameters{T}((300,300),(150.,150.))
S=generate_parameters(object_params,tau)


data_params=DatasetParameters((256,512), 64, 2,8, (T.(128.),T.(128.)))

indices=get_indices_table(data_params)
polar_params=set_default_polarisation_coefficients(T,indices)

field_params=FieldTransformParameters[]
for i=1:data_params.frames_total
    push!(field_params, FieldTransformParameters(ker,
                                                0.,
                                                (0.,0.),
                                                (-10.7365 , 1.39344),
                                                polar_params[i][1],
                                                polar_params[i][2]))
end

field_transforms=load_field_transforms(object_params,
                                       data_params,
                                       field_params)

	
psf_center=readdlm("data/PSF_centers_Airy.txt",T);
psf=readfits(Array{T,2},"data/PSF_parametered_Airy.fits");
blur=set_fft_operator(object_params,(psf[1:end÷2,:]'), psf_center[1:2])[1];


    
H = LinearDirectModel(size(S), (256,512,64),S.parameter_type,field_transforms,blur)
typeof(H*S)
    

BadPixMap=T.(rand(0.0:1e-16:1.0,data_params.size).< 0.9);

data, weight = data_simulator(BadPixMap, field_transforms, S,A=blur);

    writefits("test_results/DATA_$(tau)_$(data_params.size[1]).fits",
          ["TYPE" => "data"],
          mapslices(transpose,data,dims=[1,2]), overwrite=true)
    writefits("test_results/WEIGHT_$(tau)_$(data_params.size[1]).fits", 
          ["TYPE" => "weights"],
          mapslices(transpose,weight,dims=[1,2]), overwrite=true)

    write(S, "test_results/TRUE_$(tau)_$(data_params.size[1]).fits")
    S_convolved = PolarimetricMap("stokes", cat(blur*S.I, blur*S.Q, blur*S.U, dims=3)) 
    write(S_convolved, "test_results/TRUE_convolved_$(tau)_$(data_params.size[1]).fits")


