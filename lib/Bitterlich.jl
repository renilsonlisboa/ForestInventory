module Bitterlich

import DataFrames

export CalcBitterlich

    function CalcBitterlich(Dados, selected_num_arvores, selected_dap_column, selected_distance_column, selected_vol_column, FAB)

        FAB = Float64(Meta.parse(FAB))

        Dados[:, :Distância_Critica] = round.((50 * (Dados[!, selected_dap_column]./100))/sqrt(FAB), digits = 2)
 
        Arvores_duvidosas = filter(row -> !ismissing(row[selected_distance_column]), Dados[:, [Symbol("$selected_num_arvores"), Symbol("$selected_distance_column"), :Distância_Critica]])
        Arvores_duvidosas[:, :Observação] = ifelse.(Arvores_duvidosas[:, Symbol(selected_distance_column)] .<= Arvores_duvidosas[:, :Distância_Critica], "Distância calculada MAIOR que a medida", "Distância calculada MENOR que a medida")
        Arvores_duvidosas[:, :Conclusão] = ifelse.(Arvores_duvidosas[:, Symbol(selected_distance_column)] .<= Arvores_duvidosas[:, :Distância_Critica], "Árvore incluida na amostragem", "Árvore excluida da amostragem")

        Dados = filter(row -> ismissing(row[Symbol("$selected_distance_column")]) || row[Symbol("$selected_distance_column")] <= row[:Distância_Critica], Dados)

        m = size(Dados,1)

        AreaBasal = m * FAB

        Dados[:, :N_Arvores] = (round.(FAB./((pi.*(Dados[!, selected_dap_column].^2))/40000), digits = 4))
        
        N_Arvores = round(sum(Dados[:, :N_Arvores]), digits = 0)
 
        Volume = sum(round.(Dados[:, :N_Arvores].*Dados[:, Symbol(selected_vol_column)], digits = 4))
        
        Dados_Filtrados = Dados[:, [1, 2, 4, 5, 7]]

        Arvores_duvidosas = Arvores_duvidosas[:, [1,4,5]]

        Estimativas = DataFrames.DataFrame(Varíaveis = ["Número de Árvores (ha)", "Área Basal (m²/ha)", "Volume (m³/ha)"], Valores =[N_Arvores, AreaBasal, Volume])

        Resultados = [Estimativas, Arvores_duvidosas, Dados_Filtrados]

        return Resultados
    end

end