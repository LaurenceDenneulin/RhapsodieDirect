using Revise
using RhapsodieDirect
using DelimitedFiles
using EasyFITS
using InterpolationKernels

if prod(readdir() .!= "test_results")     
    mkdir("test_results")
end

ker=CatmullRomSpline(Float64, Flat)
tau=0.25
object_params=ObjectParameters((300,300),(150.,150.))
S=generate_parameters(object_params,tau)


data_params=DatasetParameters((256,512), 64, 2,8, (128.,128.))

indices=get_indices_table(data_params)
polar_params=set_default_polarisation_coefficients(indices)

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

	
psf_center=readdlm("data/PSF_centers_Airy.txt");
psf=readfits("data/PSF_parametered_Airy.fits");
blur=set_fft_operator(object_params,(psf[1:end÷2,:]'), psf_center[1:2])[1];


BadPixMap=Float64.(rand(0.0:1e-16:1.0,data_params.size).< 0.9);

data, weight, S_convolved=data_simulator(BadPixMap, field_transforms, blur, S);

    writefits("test_results/DATA_$(tau)_$(data_params.size[1]).fits",
          ["TYPE" => "data"],
          mapslices(transpose,data,dims=[1,2]), overwrite=true)
    writefits("test_results/WEIGHT_$(tau)_$(data_params.size[1]).fits", 
          ["TYPE" => "weights"],
          mapslices(transpose,weight,dims=[1,2]), overwrite=true)

    write(S, "test_results/TRUE_$(tau)_$(data_params.size[1]).fits")
    write(S_convolved, "test_results/TRUE_convolved_$(tau)_$(data_params.size[1]).fits")


