#!/usr/bin/env julia

using DataFrames
using PrettyTables
using CSV

in_dir = "results/OptimalDTSPN"
timelimit = 3600

function loadCSV(file)
    !isfile(file) && return nothing
    CSV.File(file) |> DataFrame
end

data = DataFrame()
for run in readdir(in_dir)
    csv_file = joinpath(in_dir, run, "data.csv")
    !isfile(csv_file) && continue
    rep = replace(run, ".txt" => "")
    type, n, id, rad, uniform = split(rep, "-")
    n, id, rad, uniform = parse.(Int, [n, id, rad, uniform])
    rad /= 100.
    csv = loadCSV(csv_file)
    filter!(cols -> cols.TIME <= timelimit, csv)
    #@show csv[end, :]
    row = csv[end, :]
    maxo = maximum(csv.OPEN)
    info = (
        Type = uniform == 1 ? "Uniform" : "Variable",    
        n = n,
        Number = id, 
        rho = rad,
        delta = rad,
        Time = row.TIME,
        LB = row.LB,
        UB = row.UB,
        GAP = row.GAP,
        var"Max|O|" = maxo,
        TreeSize = row.TREE,
    )
    push!(data, info)
    #break
end

sort!(data, [:Type, :n, :rho, :Number])

open("summary/OptimalDTSPN.html", "w") do io
    pretty_table(
        io,
        data;
        backend = Val(:html),
        tf = tf_html_minimalist
    )
end

open("summary/OptimalDTSPN.md", "w") do io
    pretty_table(
        io,
        data;
        tf = tf_markdown
    )
end

CSV.write("summary/OptimalDTSPN.csv", data)