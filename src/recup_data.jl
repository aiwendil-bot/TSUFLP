using CSV, DataFrames, DelimitedFiles

#=

get terminals : écrit dans un fichier dans un sous dossier de data (nom du sous dossier en params)
les terminaux situés dans un rectangle défini par les params (lat_inf, etc)
les terminaux sont récupérés depuis le dataset communes


get concentrators : récupère les points d'un dataset obtenu via overpass-turbo et 
les écrit dans un fichier (sous dossier de data, nom en param)
=#

function get_terminals(data::String, out::String,lat_inf::Float64,lat_sup::Float64,long_inf::Float64,long_sup::Float64)

    #conversion
    df = DataFrame(CSV.File(data))
    coordinates = df[ :,["imb_x","imb_y"]]

    e = 2.7182818284
    X = 20037508.34

    #const long4326 = (lat3857*180)/X
    coordinates.longitude = (coordinates.imb_x .* 180)  ./ X


    coordinates.latitude = coordinates.imb_y ./ (X / 180)

    exponent = (pi / 180) .* coordinates.latitude

    coordinates.latitude = map(x -> atan(e^x),exponent)
    coordinates.latitude = coordinates.latitude ./ (pi / 360)
    coordinates.latitude = coordinates.latitude .- 90

    box = coordinates[([ long_inf <= i <= long_sup for i in coordinates.longitude]) .& [ lat_inf <= i <= lat_sup for i in coordinates.latitude],:]

    if !ispath("data/$out")
        mkdir("data/$out")
    end
       open("data/$out/terminaux.txt", "w") do io
           writedlm(io, [[box.latitude[i],box.longitude[i]] for i in 1:size(box,1)])
       end

end



function get_concentrators(fileclvl1::String,fileclvl2::String,out::String)

    function string_to_mean(string::String)::Vector{Float64}

        tuple = Tuple(parse.(Float64, split(string, ',')))
        return [sum([tuple[i] for i in 2:2:length(tuple)])/(length(tuple)/2),
                sum([tuple[i] for i in 1:2:length(tuple)])/(length(tuple)/2)]
    end

    for file in [fileclvl1,fileclvl2]

        df = DataFrame(CSV.File(file))
        polygons = df[ [i == "Polygon" for i in df.geometry], ["coordinates"]]
        points = df[ [i == "Point" for i in df.geometry], ["latitude","longitude"]]
        vect_points = Vector{Vector{Float64}}(undef,size(points,1))
        for i in 1:size(points,1)
            vect_points[i] = [points[i,1],points[i,2]]
        end

        

        polygons.coordinates = map(elem->string_to_mean(elem), polygons.coordinates)

        if !ispath("data/$out")
            mkdir("data/$out")
        end

        if file == fileclvl1

            open("data/$out/clvl1.txt", "w") do io
                writedlm(io, vcat(polygons.coordinates,vect_points))
            end
        else
            open("data/$out/clvl2.txt", "w") do io
                writedlm(io, vcat(polygons.coordinates,vect_points))
            end
            
        end

    end
end
