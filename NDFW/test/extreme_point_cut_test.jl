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
        #"Braess-Example",
        #"SiouxFalls",
        "Berlin-Mitte-Center",
        #"bmc_simp_2"
            ]
            
vals_to_achieve = Dict()

vals_to_achieve["Braess-Example"] = [
    409.83,
]

vals_to_achieve["SiouxFalls"] = [
    4.67674439e6,
]

vals_to_achieve["Berlin-Mitte-Center"] = [
    #1.04961419e6
    1.10198659e6
]

ds = "Berlin-Friedrichshain"
#for ds in datasets
    tn = NDFW.load_ta_network(ds)

    _ , _, res_ta = NDFW.ta_frank_wolfe(tn)


    sres = res_ta / 100 #used for scaling objective
    cost_of_expansion = (res_ta / tn.number_of_links) #cost of expansion set to average cost of flow across arcs

    #sres = 1.0
    #cost_of_expansion = 3.543770959346052
    #output_fw = NDFW.network_design(ds, "FW", cost_of_expansion, sres)
    #output_nlmo = NDFW.network_design(ds, "NLMO-FW", cost_of_expansion, sres)

    edges2remove = [(255, 5)]
    output_ifw, x_ifw, rem_ifw = NDFW.network_design(ds, "IFW", cost_of_expansion, sres)
    output_nifw, x_nifw, rem_nifw = NDFW.network_design(ds, "NLMO-IFW", cost_of_expansion, sres)


    
    @info "output $(output_ifw)"
    @info "output $(output_nifw)"
    @info "optimal value ifw $(round(sres*output_ifw[3],digits=sig_dig))"
    @info "optimal value nifw $(round(sres*output_nifw[3],digits=sig_dig))"

    #r_ifw_added = [idx for (idx,x) in enumerate(x_ifw[713:730]) if x == 1]
    #r_nifw_added = [idx for (idx,x) in enumerate(x_nifw[713:730]) if x == 1]

    #@info "ifw added $(r_ifw_added)"
    #@info "nifw added $(r_nifw_added)"


    #@testset "Network design for dataset $(ds)" begin
        #@test res_ta ≈ vals_to_achieve[ds][1] atol=1e-2
        #@test round(sres*output_fw[3],digits=sig_dig) ≈ vals_to_achieve[ds][2] rtol=1e-4
        #@test round(sres*output_nlmo[3],digits=sig_dig) ≈ vals_to_achieve[ds][3] rtol=1e-4
    #    @test round(sres*output_nifw[3],digits=sig_dig) ≈ vals_to_achieve[ds][1] rtol=1e-4
    #end
#end

