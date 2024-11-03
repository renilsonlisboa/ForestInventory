module TCCFunctions

# Write your package code here.
using Genie, GenieFramework, DataFrames, CSV, Random, PlotlyBase
include(joinpath(@__DIR__, "app.jl")) # Amostragem Aleatória Simples
@genietools

export RunApp

    function RunApp()
        println(dirname(@__FILE__))
        Genie.loadapp(dirname(@__FILE__))
        up()
    end

end
