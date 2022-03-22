include("Libraries/ECMWFAnalysis.jl")

function CDD_HN(data_folder, out_folder, initial_year_1, initial_year_2, final_year_1, final_year_2; time_step=1)

    years_1 = initial_year_1 : time_step : final_year_1
    years_2 = initial_year_2 : time_step : final_year_2

    infilenames_1 = [string(data_folder, y, "_cold_1.grib") for y in years_2]

    infilenames_2 = [string(data_folder, y, "_cold_2.grib") for y in years_1]

    outfilenames = [string(out_folder, y, "_cold.txt") for y in years_1]

    @inbounds for i in 1 : length(outfilenames)

        println(string("Winter ", years_1[i], "-", years_2[i]))

        @time compute_CDD_HN(infilenames_1[i], infilenames_2[i], outfilenames[i])

    end

end

function MGDD_HN(data_folder, out_folder, initial_year, final_year; time_step=1)

    years = initial_year : time_step : final_year

    infilenames = [string(data_folder, y, "_hot.grib") for y in years]

    outfilenames = [string(out_folder, y, "_hot.txt") for y in years]

    @inbounds for i in 1 : length(outfilenames)

        println(string("Summer ", years[i]))

        @time compute_MGDD_HN(infilenames[i], outfilenames[i])

    end

end

function CDD_HS(data_folder, out_folder, initial_year, final_year; time_step=1)

    years = initial_year : time_step : final_year

    infilenames = [string(data_folder, y, "_cold.grib") for y in years]

    outfilenames = [string(out_folder, y, "_cold.txt") for y in years]

    @inbounds for i in 1 : length(outfilenames)

        println(string("Winter ", years[i]))

        @time compute_CDD_HS(infilenames[i], outfilenames[i])

    end

end

function MGDD_HS(data_folder, out_folder, initial_year_1, initial_year_2, final_year_1, final_year_2; time_step=1)

    years_1 = initial_year_1 : time_step : final_year_1
    years_2 = initial_year_2 : time_step : final_year_2

    infilenames_1 = [string(data_folder, y, "_hot_1.grib") for y in years_2]

    infilenames_2 = [string(data_folder, y, "_hot_2.grib") for y in years_1]

    outfilenames = [string(out_folder, y, "_hot.txt") for y in years_1]

    @inbounds for i in 1 : length(outfilenames)

        println(string("Summer ", years_1[i], "-", years_2[i]))

        @time compute_MGDD_HS(infilenames_1[i], infilenames_2[i], outfilenames[i])

    end

end

initial_year_1 = 2019
initial_year_2 = 2020

final_year_1 = 2019
final_year_2 = 2020

initial_year = 2019
final_year = 2019

println("Computing CDD")
@time CDD_HN("Data/Europe_", "Data/Europe_", initial_year_1, initial_year_2, final_year_1, final_year_2)

println("Computing MGDD")
@time MGDD_HN("Data/Europe_", "Data/Europe_", initial_year, final_year)

println("Finished!")

#@time CDD_HS(data_folder="/data/geo/ecmwf/temp_2m/Argentina/Argentina_", out_folder="GDD_data/Argentina/Argentina_", initial_year, final_year)
#@time MGDD_HS(data_folder="/data/geo/ecmwf/temp_2m/Argentina/Argentina_", out_folder="GDD_data/Argentina/Argentina_", initial_year_1, initial_year_2, final_year_1, final_year_2)


