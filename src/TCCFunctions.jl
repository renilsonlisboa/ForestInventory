module TCCFunctions

# Write your package code here.
using GenieFramework, DataFrames, CSV, Random, PlotlyBase

export RunApp

    function RunApp()
        Genie.loadapp("src/")
        up()
    end

end
