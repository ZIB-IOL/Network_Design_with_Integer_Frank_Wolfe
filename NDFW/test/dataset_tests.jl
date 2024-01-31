#This script runs the code on data sets from the Transportation Networks library and 
#compares the results to what is expected.

using DataFrames
using NDFW
using Random

seed = 0
@show seed
Random.seed!(seed)

sig_dig = 2
datasets = [
        "Braess-Example",
        "SiouxFalls",
        "Berlin-Mitte-Center",
            ]
            
vals_to_achieve = Dict()

vals_to_achieve["Braess-Example"] = [
    386.0,
    386.0,
    386.0,
    409.83,
]

vals_to_achieve["SiouxFalls"] = [
    4.2313406e6,
    4.23134116e6,
    4.23137943e6,
    4.67674439e6,
]

vals_to_achieve["Berlin-Mitte-Center"] = [
    992954.88,
    992954.7,
    992954.7,
    1.04961419e6
]


for ds in datasets
    tn = NDFW.load_ta_network(ds)

    _ , _, res_ta = NDFW.ta_frank_wolfe(tn)

    sres = res_ta / 100 #used for scaling objective
    cost_of_expansion = (res_ta / tn.number_of_links) #cost of expansion set to average cost of flow across arcs

    output_fw = NDFW.network_design(ds, "FW", cost_of_expansion, sres)
    output_nlmo = NDFW.network_design(ds, "NLMO-FW", cost_of_expansion, sres)
    output_ifw = NDFW.network_design(ds, "IFW", cost_of_expansion, sres)
    
    @testset "Network design for dataset $(ds)" begin
        @test res_ta ≈ vals_to_achieve[ds][1] atol=1e-2
        @test round(sres*output_fw[3],digits=sig_dig) ≈ vals_to_achieve[ds][2] rtol=1e-4
        @test round(sres*output_nlmo[3],digits=sig_dig) ≈ vals_to_achieve[ds][3] rtol=1e-4
        @test round(sres*output_ifw[3],digits=sig_dig) ≈ vals_to_achieve[ds][4] rtol=1e-4
    end
end

