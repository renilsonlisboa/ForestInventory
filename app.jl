# the app.jl file is the main entry point to the application. It is a bridge between the UI and the data processing logic.
using GenieFramework, DataFrames, CSV, Random, PlotlyBase
@genietools

# Define onde serão alocados os arquivos upados
const FILE_PATH = joinpath("public", "uploads")
mkpath(FILE_PATH)


# in the reactive code block, we'll implement the logic to handle user interaction
@app begin
    # we first declare reactive variables to hold the state of the interactive UI components
    # for the select menu, we need a list of station codes, and a variable to hold the selected station code
    @out method_swampling = ["Área Fixa", "Bitterlich", "Prodan", "Strand", "3P"]
    @in selected_method = ""

    # same for the metrics
    @out process_swampling = ["Amostragem Aleatória Simples", "Amostragem Estratificada", "Amostragem Sistemática", "Amostragem em Dois Estágios", "Amostragem em Conglomerados", "Amostragem Sistemática com Múltiplos Inícios Aleatórios", "Amostragem Independente", "Amostragem com Repetição Total", "Amostragem com Repetição Dupla", "Amostragem com Repetição Parcial"]
    @in selected_process = ""

    @out variables = names(CSV.read("exemplo/de.csv", DataFrame))
    @in selected_variable = ""

    @out upload_datasets = readdir("public/uploads/")
    @in selected_dataset = ""

    @in plot_area = ""
    @in inventoried_area = ""
    @in ear = ""
    @in alpha = ""
    @in selected_unit_ad = ""
    @in selected_sub_unit = ""
    @in selected_ocasiao_1 = ""
    @in selected_ocasiao_2 = ""
    @in inventoried_area_occasion_1 = ""
    @in inventoried_area_occasion_2 = ""
    @in n_potencial = ""
    @in Button_process = false
    @in selected_estrato = ""
    @in selected_subestrato = ""
    @in selected_dap_column = ""
    @in selected_vol_column = ""
    @in selected_num_arvores_column = ""
    @in fab = ""
    @in selected_distance_column = ""
    @in Distância = ""

    @out unit = ["m³", "Kg", "t", "arv"]
    @in selected_unit = ""

    # Definição do processo de importação de arquivos
    @out upfiles = readdir(FILE_PATH)
    @onchange fileuploads begin
        if !isempty(fileuploads)
            @info "Arquivo importado: " fileuploads
            filename = fileuploads["name"]
            upload_datasets = readdir("public/uploads/")
            try
                isdir(FILE_PATH) || mkpath(FILE_PATH)
                mv(fileuploads["path"], joinpath(FILE_PATH, filename), force=true)
            catch e
                @error "Error processing file: $e"
                notify(__model__, "Error processing file: $(fileuploads["name"])")
            end

            fileuploads = Dict{AbstractString,AbstractString}()
        end
        upfiles = readdir(FILE_PATH)
        upload_datasets = readdir("public/uploads/")
    end
    @event uploaded begin
        @info "uploaded"
        notify(__model__, "Arquivo Importado")
    end
    @event rejected begin
        @info "rejected"
        notify(__model__, "Falha ao Importar")
    end

    @in showUploadModal = true
    @out Table_data = DataTable()
    @out Result_Table = DataTable()
    @out Result_Table_Arvores_Duvidosas = DataTable()
    @out Resulta_Table_Filter_Data = DataTable()

    @onchange selected_dataset begin
        if selected_dataset === ""
            variables = names(CSV.read("exemplo/de.csv", DataFrame))
        else
            variables = names(CSV.read("public/uploads/$(selected_dataset)", DataFrame))
            Table_data = DataTable(CSV.read("public/uploads/$(selected_dataset)", DataFrame))
        end
    end

    # Define a visibilidade das DIVS    
    @in showUploader = true
    @in visibility_visual = false
    @in selected_unit_visibility = false
    @in selected_variable_visibility = false
    @in selected_process_visibility = false
    @in selected_dataset_visibility = false
    @in selected_unit_dropbox_visibility = false
    @in selected_subunit_dropbox_visibility = false
    @in selected_occasion_visibility = false
    @in plot_area_visibility = false
    @in inventoried_area_visibility = false
    @in inventoried_area_occasion_1_visibility = false
    @in inventoried_area_occasion_2_visibility = false
    @in ear_visibility = false
    @in alpha_visibility = false
    @in fab_visibility = false
    @in n_potencial_visibility = false
    @in selected_estrato_visibility = false
    @in selected_subestrato_visibility = false
    @in selected_dap_column_visibility = false
    @in selected_distance_column_visibility = false
    @in selected_vol_column_visibility = false
    @in selected_num_arvores_column_visibility = false
    @in Distância_visibility = false

    
    # Apresenta ou oculta os dropbox "Processo de Amostragem" e "Formato da Parcela" quando selecionado
    @onchange selected_method begin
        if selected_method === "Área Fixa" || selected_method === ""
            selected_process_visibility = true
            fab_visibility = false
            selected_dataset_visibility = false
            plot_area_visibility = false
            n_potencial_visibility = false
            selected_dap_column_visibility = false
            selected_distance_column_visibility = false
            selected_vol_column_visibility = false
            selected_num_arvores_column_visibility = false
            Distância_visibility = false
            selected_estrato_visibility = false
            selected_subestrato_visibility = false
        elseif selected_method === "Bitterlich"
            selected_unit_visibility = false
            selected_variable_visibility = false
            selected_process_visibility = false
            selected_dataset_visibility = true
            selected_unit_dropbox_visibility = false
            selected_subunit_dropbox_visibility = false
            selected_occasion_visibility = false
            plot_area_visibility = false
            inventoried_area_visibility = false
            ear_visibility = false
            alpha_visibility = false
            fab_visibility = true
            selected_dap_column_visibility = true
            selected_distance_column_visibility = true
            selected_vol_column_visibility = true
            selected_num_arvores_column_visibility = true
            Distância_visibility = false
            selected_estrato_visibility = false
            selected_subestrato_visibility = false
        elseif selected_method === "Prodan"
            selected_unit_visibility = false
            selected_variable_visibility = false
            selected_process_visibility = false
            selected_dataset_visibility = true
            selected_unit_dropbox_visibility = false
            selected_subunit_dropbox_visibility = false
            selected_occasion_visibility = false
            plot_area_visibility = false
            inventoried_area_visibility = false
            ear_visibility = false
            alpha_visibility = false
            fab_visibility = false
            selected_dap_column_visibility = true
            selected_distance_column_visibility = false
            selected_vol_column_visibility = true
            selected_num_arvores_column_visibility = false
            Distância_visibility = true
            selected_estrato_visibility = false
            selected_subestrato_visibility = false
        elseif selected_method === "Strand"
            selected_unit_visibility = false
            selected_variable_visibility = false
            selected_process_visibility = false
            selected_dataset_visibility = true
            selected_unit_dropbox_visibility = false
            selected_subunit_dropbox_visibility = false
            selected_occasion_visibility = false
            selected_dap_column_visibility = true
            selected_distance_column_visibility = false
            selected_vol_column_visibility = false
            selected_num_arvores_column_visibility = false
            plot_area_visibility = false
            inventoried_area_visibility = false
            ear_visibility = false
            alpha_visibility = false
            fab_visibility = true
            selected_estrato_visibility = false
            selected_subestrato_visibility = false
        elseif selected_method === "3P"
            selected_unit_visibility = false
            selected_variable_visibility = false
            selected_process_visibility = false
            selected_dataset_visibility = true
            selected_unit_dropbox_visibility = false
            selected_subunit_dropbox_visibility = false
            selected_occasion_visibility = false
            plot_area_visibility = true
            inventoried_area_visibility = false
            ear_visibility = false
            alpha_visibility = false
            fab_visibility = true
            selected_estrato_visibility = false
            selected_subestrato_visibility = false
        end
        selected_process, selected_unit, selected_dataset, selected_variable, plot_area, inventoried_area, ear, alpha = "", "", "", "", "", "", "", ""
    end

    # Alterna a visibilidade das DIVS no HTML conforme a seleção dos processos
    @onchange selected_process begin
        if selected_process === ""
            variables = names(CSV.read("exemplo/de.csv", DataFrame))
        elseif selected_process === "Amostragem Aleatória Simples"
            selected_unit_visibility = true
            selected_variable_visibility = true
            selected_dataset_visibility = true
            selected_occasion_visibility = false
            plot_area_visibility = true
            inventoried_area_visibility = true
            ear_visibility = true
            alpha_visibility = true
            inventoried_area_occasion_1_visibility = false
            inventoried_area_occasion_2_visibility = false
            n_potencial_visibility = false
            selected_estrato_visibility = false
            selected_subestrato_visibility = false
        elseif selected_process === "Amostragem com Repetição Dupla"
            selected_unit_visibility = true
            selected_variable_visibility = false
            selected_dataset_visibility = true
            
            selected_occasion_visibility = true
            plot_area_visibility = true
            inventoried_area_visibility = true
            ear_visibility = false
            alpha_visibility = true
        elseif selected_process === "Amostragem com Repetição Parcial"
            selected_unit_visibility = true
            selected_variable_visibility = false
            selected_dataset_visibility = true
            
            selected_occasion_visibility = true
            selected_subunit_dropbox_visibility = true
            selected_unit_dropbox_visibility = true
            plot_area_visibility = true
            inventoried_area_visibility = true
            ear_visibility = false
            alpha_visibility = true
        elseif selected_process === "Amostragem com Repetição Total"
            selected_unit_visibility = true
            selected_variable_visibility = false
            selected_dataset_visibility = true
            
            selected_occasion_visibility = true
            selected_subunit_dropbox_visibility = false
            selected_unit_dropbox_visibility = true
            plot_area_visibility = true
            inventoried_area_visibility = false
            ear_visibility = false
            alpha_visibility = true
            inventoried_area_occasion_1_visibility = true
            inventoried_area_occasion_2_visibility = true
        elseif selected_process === "Amostragem em Conglomerados"
            selected_unit_visibility = true
            selected_variable_visibility = true
            selected_dataset_visibility = true
            
            selected_occasion_visibility = false
            plot_area_visibility = true
            inventoried_area_visibility = true
            ear_visibility = true
            alpha_visibility = true
            selected_variable_visibility = false
        elseif selected_process === "Amostragem em Dois Estágios"
            selected_unit_visibility = true
            selected_variable_visibility = false
            selected_dataset_visibility = true
            
            selected_occasion_visibility = true
            plot_area_visibility = true
            inventoried_area_visibility = true
            ear_visibility = true
            alpha_visibility = true
            selected_unit_dropbox_visibility = false
            selected_subunit_dropbox_visibility = false
            n_potencial_visibility = true
            selected_estrato_visibility = false
            selected_subestrato_visibility = false
        elseif selected_process === "Amostragem Estratificada"
            selected_unit_visibility = true
            selected_variable_visibility = true
            selected_dataset_visibility = true
            
            selected_occasion_visibility = false
            plot_area_visibility = true
            inventoried_area_visibility = true
            ear_visibility = true
            alpha_visibility = true
            inventoried_area_occasion_1_visibility = false
            inventoried_area_occasion_2_visibility = false
            n_potencial_visibility = false
            selected_estrato_visibility = true
            selected_subestrato_visibility = true
        elseif selected_process === "Amostragem Independente"
            selected_unit_visibility = true
            selected_variable_visibility = false
            selected_dataset_visibility = true
            
            selected_occasion_visibility = true
            selected_subunit_dropbox_visibility = false
            selected_unit_dropbox_visibility = true
            plot_area_visibility = true
            inventoried_area_visibility = false
            ear_visibility = false
            alpha_visibility = true
            inventoried_area_occasion_1_visibility = true
            inventoried_area_occasion_2_visibility = true
        elseif selected_process === "Amostragem Sistemática com Múltiplos Inícios Aleatórios"
            selected_unit_visibility = true
            selected_variable_visibility = true
            selected_dataset_visibility = true
            selected_occasion_visibility = false
            plot_area_visibility = true
            inventoried_area_visibility = true
            ear_visibility = true
            alpha_visibility = true
            selected_estrato_visibility = false
            selected_subestrato_visibility = false
        elseif selected_process === "Amostragem Sistemática"
            selected_unit_visibility = true
            selected_variable_visibility = false
            selected_dataset_visibility = true
            selected_occasion_visibility = false
            plot_area_visibility = true
            inventoried_area_visibility = true
            ear_visibility = true
            alpha_visibility = true
            inventoried_area_occasion_1_visibility = false
            inventoried_area_occasion_2_visibility = false
            n_potencial_visibility = false
            selected_estrato_visibility = false
            selected_subestrato_visibility = false
        end
    end

    trace1 = scatter(
        x=[1, 2, 3, 4],
        y=[10, 15, 13, 17],
        mode="markers",
        name="Trace 1"
    )

    trace2 = scatter(
        x=[1, 2, 3, 4],
        y=[5, 9, 11, 12],
        mode="lines+markers",
        name="Trace 2",
        line=attr(color="red")
    )

    layout2 = PlotlyBase.Layout(
        title="A Scatter Plot with Multiple Traces",
        xaxis=attr(
            title="X Axis Label",
            showgrid=false
        ),
        yaxis=attr(
            title="Y Axis Label",
            showgrid=true,
            range=[0, 20]
        )
    )

    @out plotdata = [trace1, trace2]
    @out plotlayout = layout2

    # Define as funcionalidades do Button_process
    @in Button_process = false
    # Define as funcionalidades do Button_process
    @in Button_return = false
    @in visibility_result = false
    @in visibility_start_data = true

    @onbutton Button_process begin
        if selected_method === "Área Fixa"
            if selected_process === "Amostragem Aleatória Simples"
                Result_Table = DataTable(FixedArea.AAS(CSV.read("public/uploads/$(selected_dataset)", DataFrame), plot_area, inventoried_area, ear, alpha, selected_variable, selected_unit))
            elseif selected_process === "Amostragem com Repetição Dupla"
                Result_Table = DataTable(FixedArea.AD(CSV.read("public/uploads/$(selected_dataset)", DataFrame), plot_area, inventoried_area, selected_ocasiao_1, selected_ocasiao_2, alpha, selected_unit))
            elseif selected_process === "Amostragem com Repetição Parcial"
                Result_Table = DataTable(FixedArea.ARP(CSV.read("public/uploads/$(selected_dataset)", DataFrame), plot_area, inventoried_area, selected_unit_ad, selected_sub_unit, selected_ocasiao_1, selected_ocasiao_2, alpha, selected_unit))
            elseif selected_process === "Amostragem com Repetição Total"
                Result_Table = DataTable(FixedArea.ART(CSV.read("public/uploads/$(selected_dataset)", DataFrame), plot_area, inventoried_area_occasion_1, inventoried_area_occasion_2, selected_unit_ad, selected_ocasiao_1, selected_ocasiao_2, alpha, selected_unit))
            elseif selected_process === "Amostragem em Conglomerados"
                Result_Table = DataTable(FixedArea.CONGL(CSV.read("public/uploads/$(selected_dataset)", DataFrame), plot_area, inventoried_area, ear, alpha, selected_unit))
            elseif selected_process === "Amostragem em Dois Estágios"
                Result_Table = DataTable(FixedArea.DE(CSV.read("public/uploads/$(selected_dataset)", DataFrame), plot_area, inventoried_area, ear, alpha, n_potencial, selected_unit))
            elseif selected_process === "Amostragem Estratificada"
                Result_Table = DataTable(FixedArea.ESTRAT(CSV.read("public/uploads/$(selected_dataset)", DataFrame), plot_area, inventoried_area, ear, alpha, selected_estrato, selected_subestrato, selected_variable, selected_unit))
            elseif selected_process === "Amostragem Independente"
                Result_Table = DataTable(FixedArea.IND(CSV.read("public/uploads/$(selected_dataset)", DataFrame), plot_area, inventoried_area_occasion_1, inventoried_area_occasion_2, selected_unit_ad, selected_ocasiao_1, selected_ocasiao_2, alpha, selected_unit))
            elseif selected_process === "Amostragem Sistemática com Múltiplos Inícios Aleatórios"
                Result_Table = DataTable(FixedArea.MULTI(CSV.read("public/uploads/$(selected_dataset)", DataFrame), plot_area, inventoried_area, ear, alpha, selected_unit))
            elseif selected_process === "Amostragem Sistemática"
                Result_Table = DataTable(FixedArea.SIST(CSV.read("public/uploads/$(selected_dataset)", DataFrame), plot_area, inventoried_area, ear, alpha, selected_unit))
            else
            end
        elseif selected_method === "Bitterlich"
            Resultados = (Bitterlich.CalcBitterlich(CSV.read("public/uploads/$(selected_dataset)", DataFrame), selected_num_arvores_column, selected_dap_column, selected_distance_column, selected_vol_column, fab))
            Result_Table = DataTable(Resultados[1])
            Result_Table_Arvores_Duvidosas = DataTable(Resultados[2])
            Resulta_Table_Filter_Data = DataTable(Resultados[3])
        elseif selected_method === "Prodan"
            Result_Table = DataTable(Prodan.CalcProdan(CSV.read("public/uploads/$(selected_dataset)", DataFrame), selected_dap_column, Distância, selected_vol_column))
        elseif selected_method === "Strand"
            Result_Table = DataTable(Strand.CalcStrand(CSV.read("public/uploads/$(selected_dataset)", DataFrame), selected_dap_column, fab))
        elseif selected_method === "3P"
            println("$(selected_method)")
        else

        end
        visibility_result = true
        visibility_start_data = false
    end

    @onbutton Button_return begin
        visibility_result = false
        visibility_start_data = true
    end

end

# declare a route at / that'll render the HTML
@page("/", "app.jl.html")