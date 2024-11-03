module Strand

import DataFrames

export CalcStrand
    function CalcStrand(Dados, selected_dap_column, FAB)
        FAB = Float64(Meta.parse(FAB))

        f = 0.498
        L = 15.7
        
        G = round((sqrt(FAB)/10) * sum(Dados[:, selected_dap_column]), digits = 2)

        N = round((200*(sqrt(FAB))/L) * (sum(1.0./Dados[:, selected_dap_column])), digits = 2)

        V = round(f * (1/10) * sum(Dados[:,selected_dap_column].^2), digits = 2)

        Resultados = DataFrames.DataFrame(Varíaveis = ["Número de Árvores (ha)", "Área Basal (m²/ha)", "Volume (m³/ha)"], Valores = [N, G, V])

        return Resultados
    end
end