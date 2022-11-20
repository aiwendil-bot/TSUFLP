using CSV, DataFrames, DelimitedFiles

function get_terminals(filename::String,lat_inf::Float64,lat_sup::Float64,long_inf::Float64,long_sup::Float64)

    #conversion
    df = DataFrame(CSV.File(filename))
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

    println(describe(coordinates))

    box = coordinates[([ long_inf <= i <= long_sup for i in coordinates.longitude]) .& [ lat_inf <= i <= lat_sup for i in coordinates.latitude],:]
       open("data/terminaux.txt", "w") do io
           writedlm(io, [[box.latitude[i],box.longitude[i]] for i in 1:size(box,1)])
       end

end

#get_terminals("data/communes",47.40904,47.50549,-0.65266,-0.42212)


function get_concentrators(filename::String)

    df = DataFrame(CSV.File(filename))
    polygons = df[ [i == "Polygon" for i in df.geometry], ["coordinates"]]
    points = df[ [i == "Point" for i in df.geometry], ["latitude","longitude"]]
    vect_points = Vector{Vector{Float64}}(undef,size(points,1))
    for i in 1:size(points,1)
        vect_points[i] = [points[i,1],points[i,2]]
    end

    function string_to_mean(string::String)::Vector{Float64}

        tuple = Tuple(parse.(Float64, split(string, ',')))
        return [sum([tuple[i] for i in 2:2:length(tuple)])/(length(tuple)/2),
                sum([tuple[i] for i in 1:2:length(tuple)])/(length(tuple)/2)]
    end

    polygons.coordinates = map(elem->string_to_mean(elem), polygons.coordinates)

    open("data/clvl1.txt", "w") do io
           writedlm(io, vcat(polygons.coordinates,vect_points))
       end

end

#get_concentrators("data/clvl1.csv")
