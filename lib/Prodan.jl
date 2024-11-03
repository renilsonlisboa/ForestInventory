module Prodan

import DataFrames

export CalcProdan

    function CalcProdan(Dados, selected_dap_column, selected_distance, selected_vol_column)
        
        selected_distance = Float64(Meta.parse(selected_distance))

        R6 = (Dados[6, selected_dap_column]/100/2) + selected_distance

        Dados[:, selected_dap_column] = Dados[:, selected_dap_column]/100

        Area_basal = round((((sum(Dados[1:5, selected_dap_column].^2)) + (((Dados[6, selected_dap_column])^2)/2))/(R6^2)) * 2500, digits = 6)
        Volume = round((((sum(Dados[1:5, selected_vol_column])) + (((Dados[6, selected_vol_column]))/2))/(pi*(R6^2))) * 10000, digits = 4)
        N_Arvores = round(55000/(pi*R6^2), digits = 4)

        Estimativas = DataFrames.DataFrame(Varíaveis = ["Número de Árvores (ha)", "Área Basal (m²/ha)", "Volume (m³/ha)"], Valores =[N_Arvores, Area_basal, Volume])

        return Estimativas
    end

end