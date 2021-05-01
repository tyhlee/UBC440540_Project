using JLD
using DataFrames
using Dates
using CSV
using Printf
using RCall

function create_weekly_data(data_name,output_name)
    df = DataFrame(CSV.File(data_name))

    # """ Pre-vaccine analysis """
    tot_days = size(df)[1]
    time_interval = 7

    df_weekly = df[1:time_interval:tot_days,:]

    scale_SIR(SIR) = SIR./1e4

    transform!(df_weekly,[:S] => ByRow(scale_SIR) => :S)
    transform!(df_weekly,[:I] => ByRow(scale_SIR) => :I)
    transform!(df_weekly,[:R] => ByRow(scale_SIR) => :R)

    CSV.write(output_name,df_weekly)
end

function process_result(data_path,data_name)
    dat = load(data_path)
    key = [i for i in keys(dat)]
    dat = dat[key[1]]
    long = map(x-> dat[x,:,:],[1;2;3])
    R"saveRDS($long,$data_name)"
end

# save weekly data
covid_csv_files = readdir("data")
covid_csv_files = covid_csv_files[occursin.("covid19_bc_",covid_csv_files)]
output_names = "data/weekly/" .* covid_csv_files
covid_csv_files = "data/" .* covid_csv_files
create_weekly_data.(covid_csv_files,output_names)

# save results as rds files
results = readdir("results/covid/solution")
output_names = "results/covid/covid_rds/" .* replace.(results,".jld" => ".rds")
results = "results/covid/solution/" .* results
process_result.(results,output_names)
