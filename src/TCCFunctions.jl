module TCCFunctions

# Write your package code here.
using Genie, GenieFramework, DataFrames, CSV, Random, PlotlyBase
include(joinpath(@__DIR__, "app.jl")) # Amostragem Aleatória Simples
@genietools

export RunApp

    function RunApp()
        Genie.loadapp("")
        up()
    end

end
