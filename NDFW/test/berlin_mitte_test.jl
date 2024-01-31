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
        #"bmc_simplified",
        "Berlin-Mitte-Center_s"
            ]
            
vals_to_achieve = Dict()


vals_to_achieve["bmc_simplified"] = [
    #1.04961419e6
    1.10198659e6
]

org_res_ta = 8922.836

edges2remove = 0
#edges2remove = [(5,6), (6,7), (7,8), (8,3), (9, 10), (10, 11), (11, 12), (11, 2)]
#edges2remove = [(307, 311), (310, 313), (253, 81), (272, 256), (263, 265)]

for ds in datasets
    tn = NDFW.load_ta_network(ds)

    #_ , _, res_ta = NDFW.ta_frank_wolfe(tn)

    sres = 1.0
    #sres = res_ta / 100 #used for scaling objective
    #cost_of_expansion = (res_ta / tn.number_of_links) #cost of expansion set to average cost of flow across arcs
    cost_of_expansion = 50.0

    #output_fw = NDFW.network_design(ds, "FW", cost_of_expansion, sres)
    #output_nlmo = NDFW.network_design(ds, "NLMO-FW", cost_of_expansion, sres)
    output_ifw = NDFW.network_design(ds, "IFW", cost_of_expansion, sres, time_limit = 18000, max_fw_iter = 5000, edges2remove=edges2remove)
    output_nifw = NDFW.network_design(ds, "NLMO-IFW", cost_of_expansion, sres, time_limit = 18000, max_fw_iter = 5000, edges2remove=edges2remove)
    
    @info "output $(output_ifw)"
    @info "output $(output_nifw)"
end

