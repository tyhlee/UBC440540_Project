using DataFrames, CSV, Dates
using JLD

# import data
df_covid_raw = DataFrame(CSV.File("data/covid19-canada.csv"))

# look at BC only
df_covid_BC = select(filter(:prname => ==("British Columbia"),df_covid_raw),
                    :date,:numconf,:numactive, :numdeaths, :numrecover)

# population of BC in 2021 is about 4.9 mil
# assume s0 is 4.9 mil [estimated pop of BC in 2021]
init_pop = 4900000
compute_s0(I,R) = init_pop - (I + R)
compute_R0(D,R) = D+R
# treat deaths as recovered
transform!(df_covid_BC,[:numdeaths,:numrecover] => ByRow(compute_R0) => :numrecover)
transform!(df_covid_BC,[:numactive,:numrecover] => ByRow(compute_s0) => :S)
rename!(df_covid_BC, :numactive => :I)
rename!(df_covid_BC, :numrecover => :R)
rename!(df_covid_BC, :date => :t)
select!(df_covid_BC,:t,:S,:I,:R)

# alternative: March 19, 20201 third wave started

# select the initial date
# choose 2020-03-26; no missing data
filter!(:t => >=(Date(2020,8,10)), df_covid_BC)

# before vaccine is administered: Dec 2020
df_covid_BC_pre = filter(:t => <(Date(2021,1,11)) ,df_covid_BC)

# after vaccine
df_covid_BC_post = filter(:t => >=(Date(2021,1,11)),df_covid_BC)


# save the output
CSV.write("data/covid_19_bc_all.csv",df_covid_BC)
CSV.write("data/covid19_bc_pre_vaccine.csv",df_covid_BC_pre)
CSV.write("data/covid19_bc_post_vaccine.csv",df_covid_BC_post)

# save as JLD
# save("data/covid_bc.jld","pre",df_covid_BC_pre,"post",df_covid_BC_post)

# before vaccine is administered: Dec 2020
df_covid_BC_second_wave = filter(:t => <(Date(2021,1,19)) ,df_covid_BC)

# after vaccine
df_covid_BC_third_wave = filter(:t => >=(Date(2021,3,19)),df_covid_BC)


# save the output
CSV.write("data/covid19_bc_second_wave.csv",df_covid_BC_second_wave)
CSV.write("data/covid19_bc_third_wave.csv",df_covid_BC_third_wave)
