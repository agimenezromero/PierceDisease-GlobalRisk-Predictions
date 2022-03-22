using GRIB
using DataFrames
using Dates
using Feather
using DelimitedFiles

#There is a bug in GRIB.jl because of the new resolution format. Date is given in only days resolution 
#so it overlaps... I fix it here so cutre xd

function grib_to_feather(infilename, outfilename)
    
    df = DataFrame()
    
    GribFile(infilename) do f
        
        m = Message(f)
        
        year = parse(Int, string(m["date"])[1:4])
        month = parse(Int, string(m["date"])[5:6])
        day = parse(Int, string(m["date"])[7:8])
        
        date = DateTime(Int64(year), Int64(month), Int64(day))
        
        T = vec(m["values"])
        
        df[!, Symbol("$date")] = T
        
        for message in f

            #date = DateTime(message["date"]) #use it when bug is fixed
            
            date += Dates.Hour(1)
            
            T = vec(data(message)[3]) #Flattened T array
             
            df[!, Symbol("$date")] = T
            
        end

    end
    
    Feather.write(outfilename, df)
    
end


function compute_MGDD_HN(infilename, outfilename ; threshold_1=285.15, threshold_2=291.15, threshold_3=301.15, threshold_4=305.15, last_T=308.15)
   
    GribFile(infilename) do f
        
        m = Message(f)
        
        GDD = zeros(size(vec(m["values"])))
        
        for message in f
            
            T = vec(message["values"])
            
            for i in eachindex(T)
               
                if T[i] > threshold_1
                   
                    if T[i] < threshold_2 #First slope
                        
                        GDD[i] += 2/3 * (T[i] - threshold_1)

                    elseif (T[i] >= threshold_2) && (T[i] < threshold_3) #Second slope

                        GDD[i] +=  2/3 * (threshold_2 - threshold_1) + (T[i] - threshold_2)

                    elseif (T[i] >= threshold_3) && (T[i] < threshold_4) #Third slope


                        GDD[i] += 2/3 * (threshold_2 - threshold_1) + (threshold_3 - threshold_2) - 1.25*(T[i] - threshold_3)

                    elseif (T[i] >= threshold_4) && (T[i] < last_T) #Last slope
                        
                        GDD[i] += 2/3 * (threshold_2 - threshold_1) + (threshold_3 - threshold_2) - 1.25*(threshold_4 - threshold_3) -3*(T[i] - threshold_4)
                    
                    end
                    
                end
                
            end
            
        end

        lats = m["latitudes"]
        lons = m["longitudes"]

        shape = (size(m["values"])[2], size(m["values"])[1])
        
        f = open(outfilename, "w")

        println(f, "#GDD\tLongitudes\tLatitudes\tShape=$shape")

        for i in 1 : length(GDD)

            println(f, GDD[i]/24.0, "\t", round(lons[i], digits=4), "\t", round(lats[i], digits=4))

        end

        close(f)

    end
    
end

function compute_CDD_HS(infilename, outfilename; T_base=279.15)
   
    GribFile(infilename) do f
        
        m = Message(f)
        
        GDD = zeros(size(vec(m["values"])))
        
        for message in f
            
            T = vec(message["values"])
            
            for i in eachindex(T)
               
                if T[i] < T_base
                    
                    GDD[i] += (T_base - T[i])
                    
                end
                
            end
            
        end
        
        lats = m["latitudes"]
        lons = m["longitudes"]

        shape = (size(m["values"])[2], size(m["values"])[1])
        
        f = open(outfilename, "w")

        println(f, "#GDD\tLongitudes\tLatitudes\tshape=$shape")

        for i in 1 : length(GDD)

            println(f, GDD[i]/24.0, "\t", round(lons[i], digits=4), "\t", round(lats[i], digits=4))

        end

        close(f)
        
        return GDD

    end
    
end

function compute_CDD_HN(infilename_1, infilename_2, outfilename; T_base=279.15)
    
    #Open first file
    GribFile(infilename_1) do f
        
        m = Message(f)
        
        date = DateTime(m["date"])

        T = vec(m["values"])
        
        global GDD_1 = zeros(length(T))
        
        for message in f
            
            T = vec(message["values"])
            
            for i in eachindex(T)
               
                if T[i] < T_base
                    
                    GDD_1[i] += (T_base - T[i])
                    
                end
                
            end
            
        end

    end
    
    #Open second file
    GribFile(infilename_2) do f
        
        m = Message(f)
        
        date = DateTime(m["date"])

        T = vec(m["values"])
        
        global GDD_2 = zeros(length(T))
        
        for message in f
            
            T = vec(message["values"])
            
            for i in eachindex(T)
               
                if T[i] < T_base
                    
                    GDD_2[i] += (T_base - T[i])
                    
                end
                
            end
            
        end
        
        global lats = m["latitudes"]
        global lons = m["longitudes"]
        
        global shape = (size(m["values"])[2], size(m["values"])[1])

    end
    
    GDD = GDD_1 .+ GDD_2
    
    f = open(outfilename, "w")

    println(f, "#GDD\tLongitudes\tLatitudes\tshape=$shape")

    for i in 1 : length(GDD)

        println(f, GDD[i]/24.0, "\t", round(lons[i], digits=4), "\t", round(lats[i], digits=4))

    end
    
    close(f)

end

function compute_MGDD_HS(infilename_1, infilename_2, outfilename ; threshold_1=285.15, threshold_2=291.15, threshold_3=301.15, threshold_4=305.15, last_T=308.15)
   
    #Open first file
    GribFile(infilename_1) do f
        
        m = Message(f)
        
        T = vec(m["values"])
        
        global GDD_1 = zeros(length(T))
        
        for message in f
            
            T = vec(message["values"])
            
            for i in eachindex(T)
               
                if T[i] > threshold_1
                   
                    if T[i] < threshold_2 #First slope
                        
                        GDD_1[i] += 2/3 * (T[i] - threshold_1)

                    elseif (T[i] >= threshold_2) && (T[i] < threshold_3) #Second slope

                        GDD_1[i] +=  2/3 * (threshold_2 - threshold_1) + (T[i] - threshold_2)

                    elseif (T[i] >= threshold_3) && (T[i] < threshold_4) #Third slope


                        GDD_1[i] += 2/3 * (threshold_2 - threshold_1) + (threshold_3 - threshold_2) - 1.25*(T[i] - threshold_3)

                    elseif (T[i] >= threshold_4) && (T[i] < last_T) #Last slope
                        
                        GDD_1[i] += 2/3 * (threshold_2 - threshold_1) + (threshold_3 - threshold_2) - 1.25*(threshold_4 - threshold_3) -3*(T[i] - threshold_4)
                    
                    end
                    
                end
                
            end
            
        end

    end
    
    #Open second file
    GribFile(infilename_2) do f
        
        m = Message(f)
        
        T = vec(m["values"])
        
        global GDD_2 = zeros(length(T))
        
        for message in f
            
            T = vec(message["values"])
            
            for i in eachindex(T)
               
                if T[i] > threshold_1
                   
                    if T[i] < threshold_2 #First slope
                        
                        GDD_2[i] += 2/3 * (T[i] - threshold_1)

                    elseif (T[i] >= threshold_2) && (T[i] < threshold_3) #Second slope

                        GDD_2[i] +=  2/3 * (threshold_2 - threshold_1) + (T[i] - threshold_2)

                    elseif (T[i] >= threshold_3) && (T[i] < threshold_4) #Third slope


                        GDD_2[i] += 2/3 * (threshold_2 - threshold_1) + (threshold_3 - threshold_2) - 1.25*(T[i] - threshold_3)

                    elseif (T[i] >= threshold_4) && (T[i] < last_T) #Last slope
                        
                        GDD_2[i] += 2/3 * (threshold_2 - threshold_1) + (threshold_3 - threshold_2) - 1.25*(threshold_4 - threshold_3) -3*(T[i] - threshold_4)
                    
                    end
                    
                end
                
            end
            
        end

        global lats = m["latitudes"]
        global lons = m["longitudes"]

        global shape = (size(m["values"])[2], size(m["values"])[1])

    end
    
    GDD = GDD_1 .+ GDD_2
    
    f = open(outfilename, "w")

	println(f, "#GDD\tLongitudes\tLatitudes\tShape=$shape")

	for i in 1 : length(GDD)

	    println(f, GDD[i]/24.0, "\t", round(lons[i], digits=4), "\t", round(lats[i], digits=4))

	end

	close(f)
    
end     

function compute_avg_GDD(infilenames, outfilename)

    N = length(infilenames)

    data, header = readdlm(infilenames[1], '\t', Any, header=true)

    GDD = data[:, 1]

    shape = header[end][7:end]

    lons = data[:, 2]

    lats = data[:, 3]

    for filename in infilenames

        data, header = readdlm(filename, '\t', Any, header=true)

        GDD += data[:, 1]

    end

    GDD = GDD / N

    f = open(outfilename, "w")

    println(f, "#GDD\tlons\tlats\tshape=$shape")

    for i in 1 : length(GDD)

        println(f, GDD[i], "\t", lons[i], "\t", lats[i])

    end
    
    close(f)
    
end

function compute_min_avg_T_min(infilename, outfilename)

    GribFile(infilename) do f

        skip(f, 1)

        m = Message(f)

        year = parse(Int, string(m["date"])[1:4])
        month = parse(Int, string(m["date"])[5:6])
        day = parse(Int, string(m["date"])[7:8])

        T_min_day = vec(m["values"]) #Minimum temperature in each pixels per day
        
        pixels = length(T_min_day)
        
        Avg_T_min_month = zeros(pixels)
        
        Min_T_min_year = ones(pixels) .* 9999.0
        
        Avg_min_T_min_year = zeros(pixels)
        
        global day_count = 1.0
        global year_count = 1.0 

        for m in f

            year_2 = parse(Int, string(m["date"])[1:4])
            month_2 = parse(Int, string(m["date"])[5:6])
            day_2 = parse(Int, string(m["date"])[7:8])
            
            #For each hour change update T_min_day
            if day_2 == day
            
                T_h = vec(m["values"])
            
                for i in 1 : pixels

                    #If new hour T for the current day and pixel is lower, update T_min day
                    if T_min_day[i] > T_h[i]

                        T_min_day[i] = T_h[i]

                    end

                end
           
            elseif day_2 != day && month_2 == month

                Avg_T_min_month += T_min_day #Add min T for each pixel
                
                T_min_day = vec(m["values"]) #Update T min for new day
                
                day = day_2
                
                day_count += 1.0
                
            elseif month_2 != month && year_2 == year
                
                Avg_T_min_month += T_min_day
                
                Avg_T_min_month = Avg_T_min_month / day_count
                
                for i in 1 : pixels
                   
                    if Min_T_min_year[i] > Avg_T_min_month[i]
                        
                        Min_T_min_year[i] = Avg_T_min_month[i] #Temp min media mes mas frio
                        
                    end
                    
                end

                println("Month:", month, " Day_count:", day_count)

                month = month_2
                day = day_2
                
                day_count = 1.0   
                
                Avg_T_min_month = zeros(pixels) 
                
                T_min_day = vec(m["values"]) #Update T min for new day
                
            elseif year_2 != year #Compute Min Avg T min of the coldest year
                
                Avg_T_min_month += T_min_day
                
                Avg_T_min_month = Avg_T_min_month / day_count

                for i in 1 : pixels
                   
                    if Min_T_min_year[i] > Avg_T_min_month[i]
                        
                        Min_T_min_year[i] = Avg_T_min_month[i] #Temp min media mes mas frio
                        
                    end
                    
                end

                Avg_min_T_min_year += Min_T_min_year #Add min T for each pixel
                
                println("Year:", year)
                
                year = year_2
                month = month_2
                day = day_2
                
                year_count += 1.0
                
                day_count = 1.0 
                
                Avg_T_min_month = zeros(pixels)  
                
                T_min_day = vec(m["values"])

                Min_T_min_year = ones(pixels) .* 9999.0

            end

        end

        
        Avg_T_min_month += T_min_day

        Avg_T_min_month = Avg_T_min_month / day_count
                
        for i in 1 : pixels
                   
            if Min_T_min_year[i] > Avg_T_min_month[i]
                        
                Min_T_min_year[i] = Avg_T_min_month[i] #Temp min media mes mas frio
                        
            end
                    
        end
        
        println("Month:", month, " Day_count:", day_count)
        
        Avg_min_T_min_year += Min_T_min_year #Add min T for each pixel
                
        println("Year:", year, " Year_count:", year_count)
        
        Avg_min_T_min_year = Avg_min_T_min_year / year_count

        lats = m["latitudes"]
        lons = m["longitudes"]

        shape = (size(m["values"])[2], size(m["values"])[1])
        
        #Write file
        f = open(outfilename, "w")

        println(f, "#Tmin\tlons\tlats\tshape=$shape")

        for i in 1 : pixels

            println(f, Avg_min_T_min_year[i], "\t", round(lons[i], digits=4), "\t", round(lats[i], digits=4))

        end
            
        close(f)

    end

end

function compute_avg_T(infilename, outfilename)
    
    GribFile(infilename) do f
       
        msg = Message(f)
        
        T = vec(msg["values"])

        T_squared = T.^2
        
        count = 1
        
        for m in f
            
            count += 1
            
            T += vec(m["values"])

            T_squared += vec(m["values"]).^2
            
        end
        
        T_avg = T / count

        var = @. (T_squared / count) - (T_avg^2)

        s = @. sqrt(var/count)
        
        lats = msg["latitudes"]
        lons = msg["longitudes"]
        
        shape = (size(msg["values"])[2], size(msg["values"])[1])
        
        f = open(outfilename, "w")
        
        println(f, "#Avg Temperature\tσ\ts\tLon\tLat\tShape=$shape")
        
        for i in 1 : length(T_avg)
            
            println(f, T_avg[i], "\t", var[i], "\t", s[i], "\t", lons[i], "\t", lats[i])
            
        end
        
        close(f)
        
    end
    
end
