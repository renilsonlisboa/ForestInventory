module FixedArea

using DataFrames
using Statistics
using Distributions
using CSV

export AAS, AD, ARP, ART, CONGL, DE, ESTRAT, IND, MULTI, SIST

    function AAS(Dados, AreaParc, Area, EAR, α, selected_column, unit)

        Area = Float64(Meta.parse(Area))
        AreaParc = Float64(Meta.parse(AreaParc))
        α = Float64(Meta.parse(α))
        EAR = Float64(Meta.parse(EAR))

        N = (Area*10000)/AreaParc

        Conversor = 10000/AreaParc

        Unidades = Dados[!,1]
        Volume = (Conversor.*(Dados[!,selected_column]))

        
        AAS = DataFrames.DataFrame(Unidades = Unidades, Volume= Volume)
        
        Media = mean(Volume) #Média
        NumUni = (length(Unidades)) #Número de unidades
        Variancia = var(Volume) #Variância 
        DesvPad = std(Volume) #Desvio padrão
        FatorCorr = (1-(length(Unidades)/N)) #Fator de correção
        LE = (0.1*Media) #Limite de erro da amostragem requerido
        t = quantile(TDist(length(Unidades)-1),1-α/2) #Valor de t

        if (1-(NumUni/N)) ≥ 0.98 #f maior ou igual a 0,98 população infinita
            População = "A população avaliada é considerada infinita"   
        elseif (1-(NumUni/N)) < 0.98 #f menor que 0,98 população finita
            População = "A população avaliada é considerada finita"    
        end     
        
        Tamanho_da_amostra =   if (1-(NumUni/N)) ≥ 0.98 #f maior ou igual a 0,98 população infinita
            #População infinita. O tamanho da amostra é calculado pela seguinte equação:
            Infinita=(((t)^2)*Variancia)/(((0.1*Media))^2) 
            round(Infinita)
        elseif (1-(NumUni/N)) < 0.98 #f menor que 0,98 população finita
            #População finita. O tamanho da amostra é calculado pela seguinte equação:
            Finita=(N*((t)^2)*Variancia)/((N*(((0.1*Media))^2))+(((t)^2)*Variancia))
            round(Finita)
        end 
    
        VarMed = (Variancia/NumUni)*(1-(NumUni/N)) #Variância média
        
        ErroPad = (DesvPad/sqrt(NumUni))*sqrt((1-(NumUni/N))) #Erro padrão

        ErroPadRel = (ErroPad/mean(Volume))*100 #Erro padrão relativo

        CV = (DesvPad/Media)*100 #Coeficiente de variação

        VarMedRel = (((((sqrt(Variancia)/Media)*100))^2)/(NumUni))*(1-(NumUni/N)) #Variância média relativa
        
        #Erro de amostragem
        ErroAmostAbs = ((t*(DesvPad))/sqrt(NumUni))*sqrt((1-(NumUni/N))) #Absoluto
        ErroAmostRel = ErroAmostAbs/Media*100 #Relativo
            
        #Limite do intervalo de confiança para média 
        LII = (mean(Volume)-(t*(sqrt(var(Volume))/sqrt(length(Unidades)))*sqrt((1-(length(Unidades)/N))))) #Inferior
        LIS = (mean(Volume)+(t*(sqrt(var(Volume))/sqrt(length(Unidades)))*sqrt((1-(length(Unidades)/N))))) #Superior
        ValTotal =  ((N*mean(Volume))/Conversor) #Total da população
        
        #Limite do intervalo de confiança para o total   
        LIItotal = ((N*Media)-N*(t*(DesvPad/sqrt(NumUni))*sqrt((1-(NumUni/N)))))/Conversor #Inferior
        LIStotal = ((N*Media)+N*(t*(DesvPad/sqrt(NumUni))*sqrt((1-(NumUni/N))))) #Inferior
        ConfMin = Media-(t*(DesvPad/sqrt(NumUni))*sqrt((1-(NumUni/N)))) #Estimativa mínima de confiança
        
        #Tabela com os resultados
        if ErroAmostRel > EAR
            Observação = "Diante do exposto, conclui-se que os resultados obtidos na amostragem não satisfazem as exigências deprecisão estabelecidas para o inventário, ou seja, um erro de amostragem máximo de ±$(EAR)% da média  para confiabilidade designada. \n\nO erro estimado foi maior que o limite fixado, sendo recomendado incluir mais unidades amostrais no inventário."
        else
            Observação  = "Diante do exposto, conclui-se que os resultados obtidos na amostragem satisfazem as exigências de precisão estabelecidas para o inventário, ou seja, um erro de amostragem máximo de ±$(EAR)% da média para confiabilidade designada. \n\nO erro estimado foi menor que o limite fixado, assim as unidades amostrais são suficientes para o inventário."
        end

        Resultados = DataFrames.DataFrame(Variáveis=["Média ($(unit)/ha)", "Limite inferior do intervalo de confiança para média ($(unit)/ha)", 
        "Limite superior do intervalo de confiança para média ($(unit)/ha)", "Total da população ($(unit))", "Limite inferior do intervalo de confiança para o total ($(unit))", 
        "Limite superior do intervalo de confiança para o total ($(unit))", "Erro padrão relativo (%)", "Área da população (ha)", "Erro da amostragem absoluto ($(unit)/ha)", "Erro padrão ($(unit)/ha)", "Desvio padrão ($(unit)/ha)", 
        "Variância ($(unit)/ha)²", "Variância da média ($(unit)/ha)²", "Variância da média relativa (%)", "Coeficiente de variação (%)", "Limite de erro da amostragem requerido", "Estimativa mínima de confiança ($(unit)/ha)",
        "Fator de correção", "População", "Número total de unidades amostrais da população", 
        "Nível de significância (α)"], Valores= round.([Media, LII, LIS, ValTotal, LIItotal, LIStotal, ErroPadRel, Area, ErroAmostAbs, ErroPad, DesvPad, Variancia, VarMed, VarMedRel, CV, LE, EAR, FatorCorr, Tamanho_da_amostra, N, α], digits = 2))
     
        return Resultados
    end

    function AD(Dados, AreaParc, Area, N, selected_unit_ad, Ocasiao_1, Ocasiao_2, α, unit)
        Area = Float64(Meta.parse(Area))
        AreaParc = Float64(Meta.parse(AreaParc))
        α = Float64(Meta.parse(α))
        N = Int64(Meta.parse(N))

        Conversor=1/AreaParc
        ###Primeira ocasião####
        Unidades = Dados[:, selected_unit_ad]
        Ocasião_1 = (Conversor.*(Dados[:,Ocasiao_1]))
        Ocasião_2 = (Conversor.*(Dados[:,Ocasiao_2]))
        AD = DataFrame(Unidades = Unidades, Ocasião_1 = Ocasião_1, Ocasião_2 = Ocasião_2)

        #Média de unidades temporárias
        j=AD[!, [:Ocasião_1]]
        g=Matrix(j)
        matriz=transpose(g)
        
        global Xu = 0
        
        for i in 1:(length(Unidades)-(length(unique(Ocasião_2))))
            global Xu += matriz[i]/(length(Unidades)-(length(unique(Ocasião_2))))
        end 

        #Média de unidades permanentes
        global Xm = 0
        
        for i in ((length(Unidades)-(length(unique(Ocasião_2))))+1):length(Unidades)
            global Xm += matriz[i]/(length(unique(Ocasião_2)))
        end 

        global a = 0
        global Sxu² = 0
        
        for i in 1:(length(Unidades)-(length(unique(Ocasião_2))))
            global a += ((matriz[i].-Xu).^2)
            global Sxu² = a/((length(Unidades)-(length(unique(Ocasião_2))))-1)
        end
        
        #Variância das unidades permanentes
        global b = 0
        global Sxm² = 0
        
        for i in ((length(Unidades)-(length(unique(Ocasião_2))))+1):length(Unidades)
            global b += ((matriz[i].-Xm).^2)
            global Sxm² = b/((length(unique(Ocasião_2)))-1)
        end

        Sxu=sqrt(Sxu²) #Desvio padrão das unidades temporárias
        Sxm=sqrt(Sxm²) #Desvio padrão das unidades permanentes
        Sx1=sqrt(sum((Ocasião_1.-((((length(Unidades)-(length(unique(Ocasião_2))))/(length(Unidades))*Xu)+
        ((length(unique(Ocasião_2)))/length(Unidades))*Xm))).^2)/
        ((length(Unidades))-1)) #Desvio padrão das unidades totais 

        Primeira_ocasião = DataFrame(Variáveis=["Média das unidades amostrais totais ($(unit)/ha)", "Média das unidades amostrais temporárias ($(unit)/ha)", 
        "Média das unidades amostrais permanentes ($(unit)/ha)", "Limite inferior do intervalo de confiança para média ($(unit)/ha)", "Limite superior do intervalo de confiança para média ($(unit)/ha)", 
        "Total da população ($(unit))", "Limite inferior do intervalo de confiança para o total ($(unit))", "Limite superior do intervalo de confiança para o total ($(unit))", 
        "Área da população (ha)", "Erro da amostragem relativo (%)", "Erro da amostragem absoluto ($(unit)/ha)", "Erro padrão ($(unit)/ha)", 
        "Desvio padrão das unidades amostrais totais ($(unit)/ha)", "Desvio padrão das unidades amostrais temporárias ($(unit)/ha)", "Desvio padrão das unidades amostrais permanentes ($(unit)/ha)", 
        "Variância das unidades amostrais totais ($(unit)/ha)²", "Variância das unidades amostrais temporárias ($(unit)/ha)²", "Variância das unidades amostrais permanentes ($(unit)/ha)²", 
        "Variância da média ($(unit)/ha)²", "Limite do erro de amostragem requerido", "Tamanho da amostra", "Número total de unidades amostradas", 
        "Unidades temporárias", "Unidades permanentes", "Proporção ótima da subamostra temporária e substituída na segunda ocasião", 
        "Proporção ótima da subamostra permanente e remedida na segunda ocasião", "Nível de de significância (α)"], Valores=[(((length(Unidades)-(length(unique(Ocasião_2))))/(length(Unidades))*Xu)+
        ((length(unique(Ocasião_2)))/length(Unidades))*Xm), Xu, Xm, ((((length(Unidades)-(length(unique(Ocasião_2))))/(length(Unidades))*Xu)+
        ((length(unique(Ocasião_2)))/length(Unidades))*Xm)-((quantile(TDist((length(Unidades))-1),1-α/2))*
        (sqrt(((sum((Ocasião_1.-((((length(Unidades)-(length(unique(Ocasião_2))))/(length(Unidades))*Xu)+
        ((length(unique(Ocasião_2)))/length(Unidades))*Xm))).^2)/((length(Unidades))-1))/
        (length(Unidades)))*(1-((length(Unidades))/N)))))), ((((length(Unidades)-(length(unique(Ocasião_2))))/(length(Unidades))*Xu)+
        ((length(unique(Ocasião_2)))/length(Unidades))*Xm)+((quantile(TDist((length(Unidades))-1),1-α/2))*
        (sqrt(((sum((Ocasião_1.-((((length(Unidades)-(length(unique(Ocasião_2))))/(length(Unidades))*Xu)+
        ((length(unique(Ocasião_2)))/length(Unidades))*Xm))).^2)/((length(Unidades))-1))/(length(Unidades)))*(1-((length(Unidades))/N)))))), 
        (N*(((length(Unidades)-(length(unique(Ocasião_2))))/(length(Unidades))*Xu)+((length(unique(Ocasião_2)))/length(Unidades))*Xm)/Conversor), 
        (((N*(((length(Unidades)-(length(unique(Ocasião_2))))/(length(Unidades))*Xu)+((length(unique(Ocasião_2)))/length(Unidades))*Xm))-N*((quantile(TDist((length(Unidades))-1),1-α/2))*
        (sqrt(((sum((Ocasião_1.-((((length(Unidades)-(length(unique(Ocasião_2))))/(length(Unidades))*Xu)+((length(unique(Ocasião_2)))/length(Unidades))*Xm))).^2)/((length(Unidades))-1))/
        (length(Unidades)))*(1-((length(Unidades))/N))))))/Conversor), (((N*(((length(Unidades)-(length(unique(Ocasião_2))))/(length(Unidades))*Xu)+
        ((length(unique(Ocasião_2)))/length(Unidades))*Xm))+N*((quantile(TDist((length(Unidades))-1),1-α/2))*
        (sqrt(((sum((Ocasião_1.-((((length(Unidades)-(length(unique(Ocasião_2))))/(length(Unidades))*Xu)+
        ((length(unique(Ocasião_2)))/length(Unidades))*Xm))).^2)/((length(Unidades))-1))/(length(Unidades)))*(1-((length(Unidades))/N))))))/Conversor), 
        N, (((quantile(TDist((length(Unidades))-1),1-α/2))*(sqrt(((sum((Ocasião_1.-(((length(Unidades)-
        (length(unique(Ocasião_2))))/(length(Unidades))*Xu)+((length(unique(Ocasião_2)))/length(Unidades))*Xm)).^2)/((length(Unidades))-1))/
        (length(Unidades)))*(1-((length(Unidades))/N)))))/(((length(Unidades)-(length(unique(Ocasião_2))))/
        (length(Unidades))*Xu)+((length(unique(Ocasião_2)))/length(Unidades))*Xm))*100, ((quantile(TDist((length(Unidades))-1),1-α/2))*(sqrt(((sum((Ocasião_1.-(((length(Unidades)-
        (length(unique(Ocasião_2))))/(length(Unidades))*Xu)+((length(unique(Ocasião_2)))/length(Unidades))*Xm)).^2)/((length(Unidades))-1))/
        (length(Unidades)))*(1-((length(Unidades))/N))))), sqrt(((sum((Ocasião_1.-((((length(Unidades)-(length(unique(Ocasião_2))))/
        (length(Unidades))*Xu)+((length(unique(Ocasião_2)))/length(Unidades))*Xm))).^2)/
        ((length(Unidades))-1))/(length(Unidades)))*(1-((length(Unidades))/N))), sqrt(sum((Ocasião_1.-((((length(Unidades)-(length(unique(Ocasião_2))))/(length(Unidades))*Xu)+
        ((length(unique(Ocasião_2)))/length(Unidades))*Xm))).^2)/((length(Unidades))-1)), sqrt(Sxu²), sqrt(Sxm²), sum((Ocasião_1.-((((length(Unidades)-(length(unique(Ocasião_2))))/
        (length(Unidades))*Xu)+((length(unique(Ocasião_2)))/length(Unidades))*Xm))).^2)/((length(Unidades))-1), Sxu², Sxm², ((sum((Ocasião_1.-((((length(Unidades)-(length(unique(Ocasião_2))))/(length(Unidades))*Xu)+
        ((length(unique(Ocasião_2)))/length(Unidades))*Xm))).^2)/((length(Unidades))-1))/(length(Unidades)))*(1-((length(Unidades))/N)), 
        (0.05*(((length(Unidades)-(length(unique(Ocasião_2))))/(length(Unidades))*Xu)+((length(unique(Ocasião_2)))/length(Unidades))*Xm)), 
        ((((quantile(TDist((length(Unidades))-1),1-α/2)))^2)*(sum((Ocasião_1.-((((length(Unidades)-
        (length(unique(Ocasião_2))))/(length(Unidades))*Xu)+((length(unique(Ocasião_2)))/length(Unidades))*Xm))).^2)/((length(Unidades))-1)))/
        (((((0.05*(((length(Unidades)-(length(unique(Ocasião_2))))/(length(Unidades))*Xu)+((length(unique(Ocasião_2)))/length(Unidades))*Xm))))^2)+
        (((((quantile(TDist((length(Unidades))-1),1-α/2)))^2)*(sum((Ocasião_1.-((((length(Unidades)-(length(unique(Ocasião_2))))/(length(Unidades))*Xu)+
        ((length(unique(Ocasião_2)))/length(Unidades))*Xm))).^2)/((length(Unidades))-1)))/N)), length(Unidades), length(Unidades)-(length(unique(Ocasião_2))), 
        length(unique(Ocasião_2)), (length(Unidades)-(length(unique(Ocasião_2))))/length(Unidades), (length(unique(Ocasião_2)))/length(Unidades), 
        α]) #Tabela de resultados
            
        ###Segunda ocasião###
        j=AD[!, [:Ocasião_1]]
        g=Matrix(j)
        matriz=transpose(g)
        a = [(matriz[i].- Xm) for i in ((length(Unidades)-(length(unique(Ocasião_2))))+1):length(Unidades)]
        h = (skipmissing(Ocasião_2.-(sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2))))))
        c = sum(a.*h)
        Sxy=c/((length(unique(Ocasião_2)))-1)
        
        #Variância da regressão
        z=[(matriz[i].^2) for i in ((length(Unidades)-(length(unique(Ocasião_2))))+1):length(Unidades)]
        
        Syx²=(1/((length(unique(Ocasião_2))-2))*((sum(skipmissing(Ocasião_2.^2)))-
        (((sum(skipmissing(Ocasião_1.*Ocasião_2))^2)/sum(z))))) 

        Segunda_ocasião = DataFrame(Variáveis=["Média das unidades amostrais permanentes", "Volume médio estimado se a segunda ocasião houvesse todas unidades amostrais", 
        "Limite inferior do intervalo de confiança para média ($(unit)/ha)", "Limite superior do intervalo de confiança para média ($(unit)/ha)", "Total da população ($(unit))", 
        "Limite inferior do intervalo de confiança para o total ($(unit))", "Limite superior do intervalo de confiança para o total ($(unit))", 
        "Área da população (ha)", "Erro da amostragem relativo (%)","Erro da amostragem absoluto ($(unit)/ha)", "Erro padrão ($(unit)/ha)", "Desvio padrão ($(unit)/ha)", 
        "Variância ($(unit)/ha)²", "Variância da regressão ($(unit)/ha)²", "Variância da média ($(unit)/ha)²", "Limite do erro de amostragem requerido", 
        "Tamanho da amostra", "Número total de unidades amostradas", "Nível de significância (α)"], Valores=[sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2))), 
        (sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2))))+((Sxy/(sqrt(Sxm²*(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/
        (length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1)))))))*((sqrt(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2)))))).^2)/
        ((length(unique(Ocasião_2))-1)))))/Sxm))*(((((length(Unidades)-(length(unique(Ocasião_2))))/(length(Unidades))*Xu)+((length(unique(Ocasião_2)))/length(Unidades))*Xm))-Xm), 
        ((sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2))))+((Sxy/(sqrt(Sxm²*(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/
        (length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1)))))))*((sqrt(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2)))))).^2)/
        ((length(unique(Ocasião_2))-1)))))/Sxm))*(((((length(Unidades)-(length(unique(Ocasião_2))))/(length(Unidades))*Xu)+
        ((length(unique(Ocasião_2)))/length(Unidades))*Xm))-Xm))-((quantile(TDist((length(unique(Ocasião_2)))-1),1-α/2))*(sqrt((Syx²/(length(unique(Ocasião_2)))+
        (((sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1))))-Syx²)/(length(Unidades))))))), 
        ((sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2))))+((Sxy/(sqrt(Sxm²*(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/
        (length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1)))))))*((sqrt(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2)))))).^2)/
        ((length(unique(Ocasião_2))-1)))))/Sxm))*(((((length(Unidades)-(length(unique(Ocasião_2))))/(length(Unidades))*Xu)+((length(unique(Ocasião_2)))/length(Unidades))*Xm))-Xm))+
        ((quantile(TDist((length(unique(Ocasião_2)))-1),1-α/2))*(sqrt((Syx²/(length(unique(Ocasião_2)))+
        (((sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1))))-Syx²)/
        (length(Unidades))))))), (N*((sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2))))+((Sxy/(sqrt(Sxm²*(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/
        (length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1)))))))*((sqrt(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2)))))).^2)/
        ((length(unique(Ocasião_2))-1)))))/Sxm))*(((((length(Unidades)-(length(unique(Ocasião_2))))/
        (length(Unidades))*Xu)+((length(unique(Ocasião_2)))/length(Unidades))*Xm))-Xm))/Conversor), 
        (((N*((sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2))))+((Sxy/(sqrt(Sxm²*(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/
        (length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1)))))))*((sqrt(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/
        (length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1)))))/Sxm))*(((((length(Unidades)-
        (length(unique(Ocasião_2))))/(length(Unidades))*Xu)+((length(unique(Ocasião_2)))/length(Unidades))*Xm))-Xm)))-N*
        ((quantile(TDist((length(unique(Ocasião_2)))-1),1-α/2))*(sqrt((Syx²/(length(unique(Ocasião_2)))+
        (((sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2)))))).^2)/
        ((length(unique(Ocasião_2))-1))))-Syx²)/(length(Unidades))))))))/Conversor), (((N*((sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2))))+
        ((Sxy/(sqrt(Sxm²*(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1)))))))*
        ((sqrt(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2)))))).^2)/
        ((length(unique(Ocasião_2))-1)))))/Sxm))*(((((length(Unidades)-(length(unique(Ocasião_2))))/(length(Unidades))*Xu)+
        ((length(unique(Ocasião_2)))/length(Unidades))*Xm))-Xm)))+N*((quantile(TDist((length(unique(Ocasião_2)))-1),1-α/2))*(sqrt((Syx²/(length(unique(Ocasião_2)))+
        (((sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2)))))).^2)/
        ((length(unique(Ocasião_2))-1))))-Syx²)/(length(Unidades))))))))/Conversor), N, (((quantile(TDist((length(unique(Ocasião_2)))-1),1-α/2))*(sqrt((Syx²/(length(unique(Ocasião_2)))+
        (((sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2)))))).^2)/
        ((length(unique(Ocasião_2))-1))))-Syx²)/(length(Unidades)))))))/((sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2))))+
        ((Sxy/(sqrt(Sxm²*(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2)))))).^2)/
        ((length(unique(Ocasião_2))-1)))))))*((sqrt(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/
        (length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1)))))/Sxm))*(((((length(Unidades)-
        (length(unique(Ocasião_2))))/(length(Unidades))*Xu)+((length(unique(Ocasião_2)))/length(Unidades))*Xm))-Xm)))*100, 
        (quantile(TDist((length(unique(Ocasião_2)))-1),1-α/2))*
        (sqrt((Syx²/(length(unique(Ocasião_2)))+(((sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/
        (length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1))))-Syx²)/(length(Unidades)))))), sqrt((Syx²/(length(unique(Ocasião_2)))+(((sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/
        (length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1))))-Syx²)/(length(Unidades))))), 
        sqrt(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1)))), 
        sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1))), 
        Syx², (Syx²/(length(unique(Ocasião_2)))+(((sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/
        (length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1))))-Syx²)/(length(Unidades)))), 
        (0.05*((sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2))))+((Sxy/(sqrt(Sxm²*(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/
        (length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1)))))))*((sqrt(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2)))))).^2)/
        ((length(unique(Ocasião_2))-1)))))/Sxm))*(((((length(Unidades)-(length(unique(Ocasião_2))))/(length(Unidades))*Xu)+
        ((length(unique(Ocasião_2)))/length(Unidades))*Xm))-Xm))), ((((quantile(TDist((length(unique(Ocasião_2)))-1),1-α/2)))^2)*((length(Unidades))*(Syx²)+((length(Unidades))-
        (length(Unidades)-(length(unique(Ocasião_2)))))*((sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/
        (length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1))))-Syx²)))/((length(Unidades))*((((0.05*((sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2))))+
        ((Sxy/(sqrt(Sxm²*(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1)))))))*
        ((sqrt(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2)))))).^2)/
        ((length(unique(Ocasião_2))-1)))))/Sxm))*(((((length(Unidades)-(length(unique(Ocasião_2))))/(length(Unidades))*Xu)+
        ((length(unique(Ocasião_2)))/length(Unidades))*Xm))-Xm)))))^2)), length(unique(Ocasião_2)), α]) #Tabela de resultados

        Mudança_crescimento = DataFrame(Variáveis=["Crescimento médio ($(unit)/ha)", "Limite inferior do intervalo de confiança para média ($(unit)/ha)", 
        "Limite superior do intervalo de confiança para média ($(unit)/ha)", "Total da população ($(unit))", "Limite inferior do intervalo de confiança para o total ($(unit))", 
        "Limite superior do intervalo de confiança para o total ($(unit))", "Área da população (ha)", "Erro da amostragem relativo (%)", 
        "Erro da amostragem absoluto ($(unit)/ha)", "Erro padrão ($(unit)/ha)", "Variância média ($(unit)/ha)²"], Valores=[((sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2))))+
        ((Sxy/(sqrt(Sxm²*(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1)))))))*
        ((sqrt(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2)))))).^2)/
        ((length(unique(Ocasião_2))-1)))))/Sxm))*(((((length(Unidades)-(length(unique(Ocasião_2))))/
        (length(Unidades))*Xu)+((length(unique(Ocasião_2)))/length(Unidades))*Xm))-Xm))-((((length(Unidades)-(length(unique(Ocasião_2))))/(length(Unidades))*Xu)+
        ((length(unique(Ocasião_2)))/length(Unidades))*Xm)), (((sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2))))+
        ((Sxy/(sqrt(Sxm²*(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/
        (length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1)))))))*
        ((sqrt(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2)))))).^2)/
        ((length(unique(Ocasião_2))-1)))))/Sxm))*(((((length(Unidades)-(length(unique(Ocasião_2))))/(length(Unidades))*Xu)+
        ((length(unique(Ocasião_2)))/length(Unidades))*Xm))-Xm))-((((length(Unidades)-
        (length(unique(Ocasião_2))))/(length(Unidades))*Xu)+((length(unique(Ocasião_2)))/length(Unidades))*Xm)))-((quantile(TDist(((((length(Unidades))-1)+
        (length(unique(Ocasião_2))-1)))-1),1-α/2))*(sqrt(((((sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2)))))).^2)/
        ((length(unique(Ocasião_2))-1))))-Syx²)*((1+(Sxy/(sqrt(Sxm²*(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/
        (length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1)))))))*
        ((Sxy/(sqrt(Sxm²*(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/
        (length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1)))))))-2))/
        (((Sxy/(sqrt(Sxm²*(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/
        (length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1))))))))^2)))/(length(Unidades)))+
        (Syx²/(length(unique(Ocasião_2))))))), (((sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2))))+
        ((Sxy/(sqrt(Sxm²*(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/
        (length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1)))))))*
        ((sqrt(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2)))))).^2)/
        ((length(unique(Ocasião_2))-1)))))/Sxm))*(((((length(Unidades)-(length(unique(Ocasião_2))))/(length(Unidades))*Xu)+
        ((length(unique(Ocasião_2)))/length(Unidades))*Xm))-Xm))-((((length(Unidades)-
        (length(unique(Ocasião_2))))/(length(Unidades))*Xu)+((length(unique(Ocasião_2)))/length(Unidades))*Xm)))+((quantile(TDist(((((length(Unidades))-1)+
        (length(unique(Ocasião_2))-1)))-1),1-α/2))*(sqrt(((((sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/
        (length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1))))-Syx²)*
        ((1+(Sxy/(sqrt(Sxm²*(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/
        (length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1)))))))*
        ((Sxy/(sqrt(Sxm²*(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2)))))).^2)/
        ((length(unique(Ocasião_2))-1)))))))-2))/(((Sxy/(sqrt(Sxm²*(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/
        (length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1))))))))^2)))/(length(Unidades)))+
        (Syx²/(length(unique(Ocasião_2))))))), N*(((sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2))))+((Sxy/(sqrt(Sxm²*(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/
        (length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1)))))))*
        ((sqrt(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2)))))).^2)/
        ((length(unique(Ocasião_2))-1)))))/Sxm))*(((((length(Unidades)-(length(unique(Ocasião_2))))/
        (length(Unidades))*Xu)+((length(unique(Ocasião_2)))/length(Unidades))*Xm))-Xm))-((((length(Unidades)-(length(unique(Ocasião_2))))/(length(Unidades))*Xu)+
        ((length(unique(Ocasião_2)))/length(Unidades))*Xm))), ((N*(((sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2))))+
        ((Sxy/(sqrt(Sxm²*(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/
        (length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1)))))))*
        ((sqrt(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/
        (length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1)))))/Sxm))*(((((length(Unidades)-
        (length(unique(Ocasião_2))))/(length(Unidades))*Xu)+((length(unique(Ocasião_2)))/length(Unidades))*Xm))-Xm))-((((length(Unidades)-
        (length(unique(Ocasião_2))))/(length(Unidades))*Xu)+((length(unique(Ocasião_2)))/length(Unidades))*Xm))))-N*((quantile(TDist(((((length(Unidades))-1)+
        (length(unique(Ocasião_2))-1)))-1),1-α/2))*
        (sqrt(((((sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2)))))).^2)/
        ((length(unique(Ocasião_2))-1))))-Syx²)*((1+(Sxy/(sqrt(Sxm²*(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/
        (length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1)))))))*
        ((Sxy/(sqrt(Sxm²*(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/
        (length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1)))))))-2))/
        (((Sxy/(sqrt(Sxm²*(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/
        (length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1))))))))^2)))/(length(Unidades)))+
        (Syx²/(length(unique(Ocasião_2)))))))), ((N*(((sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2))))+
        ((Sxy/(sqrt(Sxm²*(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/
        (length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1)))))))*
        ((sqrt(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/
        (length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1)))))/Sxm))*(((((length(Unidades)-
        (length(unique(Ocasião_2))))/(length(Unidades))*Xu)+((length(unique(Ocasião_2)))/length(Unidades))*Xm))-Xm))-((((length(Unidades)-
        (length(unique(Ocasião_2))))/(length(Unidades))*Xu)+((length(unique(Ocasião_2)))/length(Unidades))*Xm))))+N*((quantile(TDist(((((length(Unidades))-1)+
        (length(unique(Ocasião_2))-1)))-1),1-α/2))*
        (sqrt(((((sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2)))))).^2)/
        ((length(unique(Ocasião_2))-1))))-Syx²)*((1+(Sxy/(sqrt(Sxm²*(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/
        (length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1)))))))*
        ((Sxy/(sqrt(Sxm²*(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/
        (length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1)))))))-2))/
        (((Sxy/(sqrt(Sxm²*(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/
        (length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1))))))))^2)))/(length(Unidades)))+
        (Syx²/(length(unique(Ocasião_2)))))))), N, (((quantile(TDist(((((length(Unidades))-1)+(length(unique(Ocasião_2))-1)))-1),1-α/2))*
        (sqrt(((((sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2)))))).^2)/
        ((length(unique(Ocasião_2))-1))))-Syx²)*((1+(Sxy/(sqrt(Sxm²*(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/
        (length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1)))))))*
        ((Sxy/(sqrt(Sxm²*(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/
        (length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1)))))))-2))/
        (((Sxy/(sqrt(Sxm²*(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/
        (length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1))))))))^2)))/(length(Unidades)))+
        (Syx²/(length(unique(Ocasião_2)))))))/(((sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2))))+
        ((Sxy/(sqrt(Sxm²*(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/
        (length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1)))))))*
        ((sqrt(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2)))))).^2)/
        ((length(unique(Ocasião_2))-1)))))/Sxm))*(((((length(Unidades)-(length(unique(Ocasião_2))))/(length(Unidades))*Xu)+
        ((length(unique(Ocasião_2)))/length(Unidades))*Xm))-Xm))-((((length(Unidades)-
        (length(unique(Ocasião_2))))/(length(Unidades))*Xu)+((length(unique(Ocasião_2)))/length(Unidades))*Xm))))*100, 
        (quantile(TDist(((((length(Unidades))-1)+(length(unique(Ocasião_2))-1)))-1),1-α/2))*
        (sqrt(((((sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2)))))).^2)/
        ((length(unique(Ocasião_2))-1))))-Syx²)*((1+(Sxy/(sqrt(Sxm²*(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/
        (length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1)))))))*
        ((Sxy/(sqrt(Sxm²*(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2)))))).^2)/
        ((length(unique(Ocasião_2))-1)))))))-2))/(((Sxy/(sqrt(Sxm²*(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/
        (length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1))))))))^2)))/(length(Unidades)))+
        (Syx²/(length(unique(Ocasião_2)))))), sqrt(((((sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2)))))).^2)/
        ((length(unique(Ocasião_2))-1))))-Syx²)*((1+(Sxy/(sqrt(Sxm²*(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/
        (length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1)))))))*((Sxy/(sqrt(Sxm²*(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/
        (length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1)))))))-2))/(((Sxy/(sqrt(Sxm²*(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/
        (length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1))))))))^2)))/(length(Unidades)))+
        (Syx²/(length(unique(Ocasião_2))))), ((((sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2)))))).^2)/
        ((length(unique(Ocasião_2))-1))))-Syx²)*((1+(Sxy/(sqrt(Sxm²*(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/
        (length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1)))))))*
        ((Sxy/(sqrt(Sxm²*(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/(length(unique(Ocasião_2)))))).^2)/
        ((length(unique(Ocasião_2))-1)))))))-2))/(((Sxy/(sqrt(Sxm²*(sum((skipmissing((Ocasião_2.-(sum(skipmissing(Ocasião_2))/
        (length(unique(Ocasião_2)))))).^2)/((length(unique(Ocasião_2))-1))))))))^2)))/(length(Unidades)))+
        (Syx²/(length(unique(Ocasião_2))))]) #Tabela de resultados

        return Resultado
    end

    function ARP(Dados, AreaParc, N, unidade, subunidade, Ocasiao_1, Ocasiao_2, α, unit)

        AreaParc = Float64(Meta.parse(AreaParc))
        α = Float64(Meta.parse(α))
        N = Int64(Meta.parse(N))

        Conversor=1/AreaParc

        ###Primeira ocasião####
        Unidades = Dados[!,unidades]
        SubAmostra = Dados[!,subunidades]
        Ocasião_1 = (Conversor.*(Dados[:,Ocasiao_1]))
        Ocasião_2 = (Conversor.*(Dados[:,Ocasiao_2]))
        ARP = DataFrame(Unidades = Unidades, SubAmostra = SubAmostra, Ocasião_1 = Ocasião_1, Ocasião_2 = Ocasião_2)

        Informações_do_inventário = DataFrame(Variáveis=["Área da população (ha)", "Número total de unidades", 
        "Número de unidades da primeira ocasião", "Número de unidades da segunda ocasião", "Número de subamostras temporárias", 
        "Número de subamostras permanentes", "Número de novas subamostra temporárias", "Proporção ótima da subamostra temporária e substituída na segunda ocasião", 
        "Proporção ótima da subamostra permanente e remedida na segunda ocasião", "Nível de significância (α)"], Valores=[N, length(Unidades), 
        length(unique(skipmissing(Ocasião_1))), length(unique(skipmissing(Ocasião_2))), (length(Unidades))-(length(unique(skipmissing(Ocasião_2)))), 
        (length(Unidades))-(length(unique(skipmissing(Ocasião_1)))), ((length(Unidades))-(length(unique(skipmissing(Ocasião_1))))), 
        ((length(Unidades))-(length(unique(skipmissing(Ocasião_2)))))/(length(unique(skipmissing(Ocasião_1)))), 
        ((length(Unidades))-(length(unique(skipmissing(Ocasião_1)))))/(length(unique(skipmissing(Ocasião_1)))), α])

        ###Primeira ocasião###  
        #Média das unidades de amostragem temporárias
        j=ARP[!, [:Ocasião_1]]
        g=Matrix(j)
        matriz=transpose(g)
        
        global Xu = 0
        
        for i in 1:((length(Unidades))-(length(unique(skipmissing(Ocasião_2)))))
            global Xu += matriz[i]/(((length(Unidades))-(length(unique(skipmissing(Ocasião_2))))))
        end
        
        #Média das unidades de amostragem permanentes
        global Xm = 0
        
        for i in (((length(Unidades))-(length(unique(skipmissing(Ocasião_2)))))+1):(length(unique(skipmissing(Ocasião_1))))
            global Xm += matriz[i]/(((length(Unidades))-(length(unique(skipmissing(Ocasião_1))))))
        end 
        
        #Variância das unidades de amostragem temporárias
        global A = 0
        global Sxu² = 0
        
        for i in 1:((length(Unidades))-(length(unique(skipmissing(Ocasião_2)))))
            global A += ((matriz[i].-Xu).^2)
            global Sxu² = A/(((length(Unidades))-(length(unique(skipmissing(Ocasião_2)))))-1)
        end
        
        #Variância das unidades de amostragem permanentes
        global C = 0
        global Sxm² = 0
        
        for i in (((length(Unidades))-(length(unique(skipmissing(Ocasião_2)))))+1):(length(unique(skipmissing(Ocasião_1))))
            global C += ((matriz[i].-Xm).^2)
            global Sxm² = C/(((length(Unidades))-(length(unique(skipmissing(Ocasião_1)))))-1)
        end

        Primeira_ocasião = DataFrame(Variáveis=["Média das unidades amostrais totais ($(unit)/ha)", "Média das unidades amostrais temporárias ($(unit)/ha)", 
        "Média das unidades amostrais permanentes ($(unit)/ha)", "Limite inferior do intervalo de confiança para média ($(unit)/ha)", 
        "Limite superior do intervalo de confiança para média ($(unit)/ha)", "Total da população ($(unit))", 
        "Limite inferior do intervalo de confiança para o total ($(unit))", "Limite superior do intervalo de confiança para o total ($(unit))", 
        "Erro da amostragem relativo (%)", "Erro da amostragem absoluto ($(unit)/ha)", "Erro padrão ($(unit)/ha)", "Variância das unidades amostrais totais ($(unit)/ha)²", 
        "Variância das unidades amostrais temporárias ($(unit)/ha)²", "Variância das unidades amostrais permanentes ($(unit)/ha)²", "Variância da média ($(unit)/ha)²"], 
        Valores=[sum(skipmissing(Ocasião_1))/(length(unique(skipmissing(Ocasião_1)))), Xu, Xm, ((sum(skipmissing(Ocasião_1))/(length(unique(skipmissing(Ocasião_1)))))-
        ((quantile(TDist((length(unique(skipmissing(Ocasião_1))))-1),1-α/2))*(sqrt(((sum((skipmissing(Ocasião_1).-(sum(skipmissing(Ocasião_1))/
        (length(unique(skipmissing(Ocasião_1)))))).^2)/(length(unique(skipmissing(Ocasião_1))))-1)/(length(unique(skipmissing(Ocasião_1)))))*
        (1-((length(unique(skipmissing(Ocasião_1))))/N)))))), ((sum(skipmissing(Ocasião_1))/(length(unique(skipmissing(Ocasião_1)))))+
        ((quantile(TDist((length(unique(skipmissing(Ocasião_1))))-1),1-α/2))*(sqrt(((sum((skipmissing(Ocasião_1).-(sum(skipmissing(Ocasião_1))/
        (length(unique(skipmissing(Ocasião_1)))))).^2)/(length(unique(skipmissing(Ocasião_1))))-1)/(length(unique(skipmissing(Ocasião_1)))))*
        (1-((length(unique(skipmissing(Ocasião_1))))/N)))))), (N*(sum(skipmissing(Ocasião_1))/(length(unique(skipmissing(Ocasião_1)))))/Conversor), 
        (((N*(sum(skipmissing(Ocasião_1))/(length(unique(skipmissing(Ocasião_1))))))-N*((quantile(TDist((length(unique(skipmissing(Ocasião_1))))-1),1-α/2))*
        (sqrt(((sum((skipmissing(Ocasião_1).-(sum(skipmissing(Ocasião_1))/(length(unique(skipmissing(Ocasião_1)))))).^2)/
        (length(unique(skipmissing(Ocasião_1))))-1)/(length(unique(skipmissing(Ocasião_1)))))*
        (1-((length(unique(skipmissing(Ocasião_1))))/N))))))/Conversor), (((N*(sum(skipmissing(Ocasião_1))/
        (length(unique(skipmissing(Ocasião_1))))))+N*((quantile(TDist((length(unique(skipmissing(Ocasião_1))))-1),1-α/2))*
        (sqrt(((sum((skipmissing(Ocasião_1).-(sum(skipmissing(Ocasião_1))/(length(unique(skipmissing(Ocasião_1)))))).^2)/
        (length(unique(skipmissing(Ocasião_1))))-1)/(length(unique(skipmissing(Ocasião_1)))))*
        (1-((length(unique(skipmissing(Ocasião_1))))/N))))))/Conversor), (((quantile(TDist((length(unique(skipmissing(Ocasião_1))))-1),1-α/2))*
        (sqrt(((sum((skipmissing(Ocasião_1).-(sum(skipmissing(Ocasião_1))/(length(unique(skipmissing(Ocasião_1)))))).^2)/
        (length(unique(skipmissing(Ocasião_1))))-1)/(length(unique(skipmissing(Ocasião_1)))))*
        (1-((length(unique(skipmissing(Ocasião_1))))/N)))))/(sum(skipmissing(Ocasião_1))/(length(unique(skipmissing(Ocasião_1))))))*100, 
        (quantile(TDist((length(unique(skipmissing(Ocasião_1))))-1),1-α/2))*(sqrt(((sum((skipmissing(Ocasião_1).-(sum(skipmissing(Ocasião_1))/
        (length(unique(skipmissing(Ocasião_1)))))).^2)/(length(unique(skipmissing(Ocasião_1))))-1)/(length(unique(skipmissing(Ocasião_1)))))*
        (1-((length(unique(skipmissing(Ocasião_1))))/N)))), sqrt(((sum((skipmissing(Ocasião_1).-(sum(skipmissing(Ocasião_1))/(length(unique(skipmissing(Ocasião_1)))))).^2)/
        (length(unique(skipmissing(Ocasião_1))))-1)/(length(unique(skipmissing(Ocasião_1)))))*
        (1-((length(unique(skipmissing(Ocasião_1))))/N))),  sum((skipmissing(Ocasião_1).-(sum(skipmissing(Ocasião_1))/(length(unique(skipmissing(Ocasião_1)))))).^2)/
        (length(unique(skipmissing(Ocasião_1))))-1, Sxu², Sxm², ((sum((skipmissing(Ocasião_1).-(sum(skipmissing(Ocasião_1))/(length(unique(skipmissing(Ocasião_1)))))).^2)/
        (length(unique(skipmissing(Ocasião_1))))-1)/(length(unique(skipmissing(Ocasião_1)))))*
        (1-((length(unique(skipmissing(Ocasião_1))))/N))]) #Tabela de resultados 
        
        #Média de unidades de amostragem permanentes
        j=ARP[!, [:Ocasião_2]]
        g=Matrix(j)
        matriz=transpose(g)
        
        global Ym = 0
        
        for i in (((length(Unidades))-(length(unique(skipmissing(Ocasião_2)))))+1):(length(unique(skipmissing(Ocasião_1))))
            global Ym += matriz[i]/(((length(Unidades))-(length(unique(skipmissing(Ocasião_1))))))
        end 
        
        #Média de novas unidades de amostragem temporárias
        j=ARP[!, [:Ocasião_2]]
        g=Matrix(j)
        matriz=transpose(g)
        
        global Yn = 0
        
        for i in ((length(unique(skipmissing(Ocasião_1))))+1):(length(Unidades))
            global Yn += matriz[i]/(((length(Unidades))-(length(unique(skipmissing(Ocasião_1))))))
        end
        
        #Variância de unidades de amostragem permanentes
        j=ARP[!, [:Ocasião_2]]
        g=Matrix(j)
        matriz=transpose(g)
        
        global B = 0
        global Sym² = 0
        
        for i in (((length(Unidades))-(length(unique(skipmissing(Ocasião_2)))))+1):(length(unique(skipmissing(Ocasião_1))))
            global B += ((matriz[i].-Ym).^2)
            global Sym² = B/(((length(Unidades))-(length(unique(skipmissing(Ocasião_1)))))-1)
        end

        #Variância de novas unidades de amostragem temporárias
        j=ARP[!, [:Ocasião_2]]
        g=Matrix(j)
        matriz=transpose(g)

        global B = 0
        global Syn² = 0
        
        for i in ((length(unique(skipmissing(Ocasião_1))))+1):(length(Unidades))
            global B += ((matriz[i].-Yn).^2)
            global Syn² = B/((((length(Unidades))-(length(unique(skipmissing(Ocasião_1))))))-1)
        end

        #Variâncias de unidades permanentes e temporárias
        j=ARP[!, [:Ocasião_1]]
        g=Matrix(j)
        matriz=transpose(g)
        D = [matriz[i].-Xm for i in (((length(Unidades))-
            (length(unique(skipmissing(Ocasião_2)))))+1):(length(unique(skipmissing(Ocasião_1))))]
        j=ARP[!, [:Ocasião_2]]
        g=Matrix(j)
        matriz=transpose(g)
        
        J = [matriz[i].-Ym for i in (((length(Unidades))-
            (length(unique(skipmissing(Ocasião_2)))))+1):(length(unique(skipmissing(Ocasião_1))))]
        sum(D.*J)/(((length(Unidades))-(length(unique(skipmissing(Ocasião_1)))))-1)
        (sum(D.*J)/(((length(Unidades))-(length(unique(skipmissing(Ocasião_1)))))-1))/
        (sqrt(Sxm²)*sqrt(Sym²)) #Coeficiente de correlação entre os volumes das unidades permanentes medidas nas duas ocasiões
        ((sum(D.*J)/(((length(Unidades))-(length(unique(skipmissing(Ocasião_1)))))-1))/
        (sqrt(Sxm²)*sqrt(Sym²)))*(sqrt(Sym²)/sqrt(Sxm²))
        
        #Constantes "a" e "b" 
        a=((((length(Unidades))-(length(unique(skipmissing(Ocasião_1)))))*(((length(Unidades))-
        (length(unique(skipmissing(Ocasião_2)))))/(length(unique(skipmissing(Ocasião_1))))))/
        ((length(unique(skipmissing(Ocasião_2))))-((((length(Unidades))-(length(unique(skipmissing(Ocasião_2)))))/
        (length(unique(skipmissing(Ocasião_1)))))*((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1)))))*((((sum(D.*J)/(((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1)))))-1))/(sqrt(Sxm²)*sqrt(Sym²))))^2))))*
        (((sum(D.*J)/(((length(Unidades))-(length(unique(skipmissing(Ocasião_1)))))-1))/(sqrt(Sxm²)*sqrt(Sym²)))*
        (sqrt(Sym²)/sqrt(Sxm²))) ##Constante a
        
        c=((length(Unidades))-(length(unique(skipmissing(Ocasião_1)))))/((length(unique(skipmissing(Ocasião_2))))-
        ((((length(Unidades))-(length(unique(skipmissing(Ocasião_2)))))/(length(unique(skipmissing(Ocasião_1)))))*
        (((length(Unidades))-(length(unique(skipmissing(Ocasião_1))))))*((((sum(D.*J)/(((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1)))))-1))/(sqrt(Sxm²)*sqrt(Sym²))))^2))) #Constante c

        Segunda_ocasião = DataFrame(Variáveis=["Média das unidades amostrais totais ($(unit)/ha)", "Média das unidades amostrais permanente ($(unit)/ha)", 
        "Média das novas unidades amostrais temporárias ($(unit)/ha)", "Média corrente ($(unit)/ha)", "Limite inferior do intervalo de confiança para média ($(unit)/ha)", 
        "Limite superior do intervalo de confiança para média ($(unit)/ha)", "Total da população ($(unit))", 
        "Limite inferior do intervalo de confiança para o total ($(unit))", "Limite superior do intervalo de confiança para o total ($(unit))", 
        "Erro da amostragem relativo (%)", "Erro da amostragem absoluto ($(unit)/ha)", "Erro padrão ($(unit)/ha)", 
        "Variância das unidades amostrais totais ($(unit)/ha)²", "Variância das novas unidades amostrais temporárias ($(unit)/ha)²", 
        "Variância das unidades amostrais permanentes ($(unit)/ha)²", "Variância da média ($(unit)/ha)²", "Coeficiente de correlação", "Constante a", 
        "Constante c"], Valores=[sum(skipmissing(Ocasião_2))/(length(unique(skipmissing(Ocasião_2)))), Yn, Ym, ((a*Xu)-(a*Xm)+(c*Ym)+(1-c)*Yn), 
        ((((a*Xu)-(a*Xm)+(c*Ym)+(1-c)*Yn))-((quantile(TDist((length(unique(skipmissing(Ocasião_2))))-1),1-α/2))*
        (sqrt(((a^2)*(sum((skipmissing(Ocasião_1).-(sum(skipmissing(Ocasião_1))/(length(unique(skipmissing(Ocasião_1)))))).^2)/
        (length(unique(skipmissing(Ocasião_1))))-1)*((1/((length(Unidades))-(length(unique(skipmissing(Ocasião_2))))))+(1/((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1))))))))+((c^2)*((sum((skipmissing(Ocasião_2).-(sum(skipmissing(Ocasião_2))/
        (length(unique(skipmissing(Ocasião_2)))))).^2)/(length(unique(skipmissing(Ocasião_2))))-1)/((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1)))))))+(((1-c)^2)*((sum((skipmissing(Ocasião_2).-(sum(skipmissing(Ocasião_2))/
        (length(unique(skipmissing(Ocasião_2)))))).^2)/(length(unique(skipmissing(Ocasião_2))))-1)/((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1)))))))-((2*a*c*((sum(D.*J)/(((length(Unidades))-(length(unique(skipmissing(Ocasião_1)))))-1))/(sqrt(Sxm²)*sqrt(Sym²))))*
        ((sqrt(Sym²)*sqrt(Sxm²))/((length(Unidades))-(length(unique(skipmissing(Ocasião_1))))))))))), 
        ((((a*Xu)-(a*Xm)+(c*Ym)+(1-c)*Yn))+((quantile(TDist((length(unique(skipmissing(Ocasião_2))))-1),1-α/2))*
        (sqrt(((a^2)*(sum((skipmissing(Ocasião_1).-(sum(skipmissing(Ocasião_1))/(length(unique(skipmissing(Ocasião_1)))))).^2)/
        (length(unique(skipmissing(Ocasião_1))))-1)*((1/((length(Unidades))-(length(unique(skipmissing(Ocasião_2))))))+(1/((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1))))))))+((c^2)*((sum((skipmissing(Ocasião_2).-(sum(skipmissing(Ocasião_2))/
        (length(unique(skipmissing(Ocasião_2)))))).^2)/(length(unique(skipmissing(Ocasião_2))))-1)/((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1)))))))+(((1-c)^2)*((sum((skipmissing(Ocasião_2).-(sum(skipmissing(Ocasião_2))/
        (length(unique(skipmissing(Ocasião_2)))))).^2)/(length(unique(skipmissing(Ocasião_2))))-1)/((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1)))))))-((2*a*c*((sum(D.*J)/(((length(Unidades))-(length(unique(skipmissing(Ocasião_1)))))-1))/
        (sqrt(Sxm²)*sqrt(Sym²))))*((sqrt(Sym²)*sqrt(Sxm²))/((length(Unidades))-(length(unique(skipmissing(Ocasião_1))))))))))), 
        (N*(((a*Xu)-(a*Xm)+(c*Ym)+(1-c)*Yn))/Conversor), (((N*(((a*Xu)-(a*Xm)+(c*Ym)+(1-c)*Yn)))-N*((quantile(TDist((length(unique(skipmissing(Ocasião_2))))-1),1-α/2))*
        (sqrt(((a^2)*(sum((skipmissing(Ocasião_1).-(sum(skipmissing(Ocasião_1))/(length(unique(skipmissing(Ocasião_1)))))).^2)/
        (length(unique(skipmissing(Ocasião_1))))-1)*((1/((length(Unidades))-(length(unique(skipmissing(Ocasião_2))))))+(1/((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1))))))))+((c^2)*((sum((skipmissing(Ocasião_2).-(sum(skipmissing(Ocasião_2))/
        (length(unique(skipmissing(Ocasião_2)))))).^2)/(length(unique(skipmissing(Ocasião_2))))-1)/((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1)))))))+(((1-c)^2)*((sum((skipmissing(Ocasião_2).-(sum(skipmissing(Ocasião_2))/
        (length(unique(skipmissing(Ocasião_2)))))).^2)/(length(unique(skipmissing(Ocasião_2))))-1)/((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1)))))))-((2*a*c*((sum(D.*J)/(((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1)))))-1))/(sqrt(Sxm²)*sqrt(Sym²))))*((sqrt(Sym²)*sqrt(Sxm²))/((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1)))))))))))/Conversor), (((N*(((a*Xu)-(a*Xm)+(c*Ym)+(1-c)*Yn)))+N*((quantile(TDist((length(unique(skipmissing(Ocasião_2))))-1),1-α/2))*
        (sqrt(((a^2)*(sum((skipmissing(Ocasião_1).-(sum(skipmissing(Ocasião_1))/(length(unique(skipmissing(Ocasião_1)))))).^2)/
        (length(unique(skipmissing(Ocasião_1))))-1)*((1/((length(Unidades))-(length(unique(skipmissing(Ocasião_2))))))+(1/((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1))))))))+((c^2)*((sum((skipmissing(Ocasião_2).-(sum(skipmissing(Ocasião_2))/
        (length(unique(skipmissing(Ocasião_2)))))).^2)/(length(unique(skipmissing(Ocasião_2))))-1)/((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1)))))))+(((1-c)^2)*((sum((skipmissing(Ocasião_2).-(sum(skipmissing(Ocasião_2))/
        (length(unique(skipmissing(Ocasião_2)))))).^2)/(length(unique(skipmissing(Ocasião_2))))-1)/((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1)))))))-((2*a*c*((sum(D.*J)/(((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1)))))-1))/(sqrt(Sxm²)*sqrt(Sym²))))*((sqrt(Sym²)*sqrt(Sxm²))/((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1)))))))))))/Conversor), (((quantile(TDist((length(unique(skipmissing(Ocasião_2))))-1),1-α/2))*(sqrt(((a^2)*(sum((skipmissing(Ocasião_1).-(sum(skipmissing(Ocasião_1))/
        (length(unique(skipmissing(Ocasião_1)))))).^2)/(length(unique(skipmissing(Ocasião_1))))-1)*((1/((length(Unidades))-
        (length(unique(skipmissing(Ocasião_2))))))+(1/((length(Unidades))-(length(unique(skipmissing(Ocasião_1))))))))+((c^2)*
        ((sum((skipmissing(Ocasião_2).-(sum(skipmissing(Ocasião_2))/(length(unique(skipmissing(Ocasião_2)))))).^2)/
        (length(unique(skipmissing(Ocasião_2))))-1)/((length(Unidades))-(length(unique(skipmissing(Ocasião_1)))))))+(((1-c)^2)*
        ((sum((skipmissing(Ocasião_2).-(sum(skipmissing(Ocasião_2))/(length(unique(skipmissing(Ocasião_2)))))).^2)/
        (length(unique(skipmissing(Ocasião_2))))-1)/((length(Unidades))-(length(unique(skipmissing(Ocasião_1)))))))-((2*a*c*((sum(D.*J)/(((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1)))))-1))/(sqrt(Sxm²)*sqrt(Sym²))))*((sqrt(Sym²)*sqrt(Sxm²))/((length(Unidades))-(length(unique(skipmissing(Ocasião_1))))))))))/
        (((a*Xu)-(a*Xm)+(c*Ym)+(1-c)*Yn)))*100, (quantile(TDist((length(unique(skipmissing(Ocasião_2))))-1),1-α/2))*
        (sqrt(((a^2)*(sum((skipmissing(Ocasião_1).-(sum(skipmissing(Ocasião_1))/(length(unique(skipmissing(Ocasião_1)))))).^2)/
        (length(unique(skipmissing(Ocasião_1))))-1)*((1/((length(Unidades))-(length(unique(skipmissing(Ocasião_2))))))+(1/((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1))))))))+((c^2)*((sum((skipmissing(Ocasião_2).-(sum(skipmissing(Ocasião_2))/(length(unique(skipmissing(Ocasião_2)))))).^2)/
        (length(unique(skipmissing(Ocasião_2))))-1)/((length(Unidades))-(length(unique(skipmissing(Ocasião_1)))))))+(((1-c)^2)*
        ((sum((skipmissing(Ocasião_2).-(sum(skipmissing(Ocasião_2))/(length(unique(skipmissing(Ocasião_2)))))).^2)/
        (length(unique(skipmissing(Ocasião_2))))-1)/((length(Unidades))-(length(unique(skipmissing(Ocasião_1)))))))-((2*a*c*((sum(D.*J)/(((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1)))))-1))/(sqrt(Sxm²)*sqrt(Sym²))))*((sqrt(Sym²)*sqrt(Sxm²))/((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1))))))))), sqrt(((a^2)*(sum((skipmissing(Ocasião_1).-(sum(skipmissing(Ocasião_1))/(length(unique(skipmissing(Ocasião_1)))))).^2)/
        (length(unique(skipmissing(Ocasião_1))))-1)*((1/((length(Unidades))-(length(unique(skipmissing(Ocasião_2))))))+(1/((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1))))))))+((c^2)*((sum((skipmissing(Ocasião_2).-(sum(skipmissing(Ocasião_2))/(length(unique(skipmissing(Ocasião_2)))))).^2)/
        (length(unique(skipmissing(Ocasião_2))))-1)/((length(Unidades))-(length(unique(skipmissing(Ocasião_1)))))))+(((1-c)^2)*
        ((sum((skipmissing(Ocasião_2).-(sum(skipmissing(Ocasião_2))/(length(unique(skipmissing(Ocasião_2)))))).^2)/
        (length(unique(skipmissing(Ocasião_2))))-1)/((length(Unidades))-(length(unique(skipmissing(Ocasião_1)))))))-((2*a*c*((sum(D.*J)/(((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1)))))-1))/(sqrt(Sxm²)*sqrt(Sym²))))*((sqrt(Sym²)*sqrt(Sxm²))/((length(Unidades))-(length(unique(skipmissing(Ocasião_1)))))))), 
        sum((skipmissing(Ocasião_2).-(sum(skipmissing(Ocasião_2))/(length(unique(skipmissing(Ocasião_2)))))).^2)/
        (length(unique(skipmissing(Ocasião_2))))-1, Syn², Sym², ((a^2)*(sum((skipmissing(Ocasião_1).-(sum(skipmissing(Ocasião_1))/(length(unique(skipmissing(Ocasião_1)))))).^2)/
        (length(unique(skipmissing(Ocasião_1))))-1)*((1/((length(Unidades))-(length(unique(skipmissing(Ocasião_2))))))+(1/((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1))))))))+((c^2)*((sum((skipmissing(Ocasião_2).-(sum(skipmissing(Ocasião_2))/(length(unique(skipmissing(Ocasião_2)))))).^2)/
        (length(unique(skipmissing(Ocasião_2))))-1)/((length(Unidades))-(length(unique(skipmissing(Ocasião_1)))))))+(((1-c)^2)*
        ((sum((skipmissing(Ocasião_2).-(sum(skipmissing(Ocasião_2))/(length(unique(skipmissing(Ocasião_2)))))).^2)/
        (length(unique(skipmissing(Ocasião_2))))-1)/((length(Unidades))-(length(unique(skipmissing(Ocasião_1)))))))-((2*a*c*((sum(D.*J)/(((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1)))))-1))/(sqrt(Sxm²)*sqrt(Sym²))))*((sqrt(Sym²)*sqrt(Sxm²))/
        ((length(Unidades))-(length(unique(skipmissing(Ocasião_1))))))), (sum(D.*J)/(((length(Unidades))-(length(unique(skipmissing(Ocasião_1)))))-1))/
        (sqrt(Sxm²)*sqrt(Sym²)), a, c]) #Tabela de resultados
            
        ###Crescimento ou mudança###
        #Estimativa direta 
        #É necessário achar os coeficientes "A" e "B" para calcular a crescimento médio 
        A=(((length(Unidades))-(length(unique(skipmissing(Ocasião_1)))))/((length(unique(skipmissing(Ocasião_2))))-
        ((((length(Unidades))-(length(unique(skipmissing(Ocasião_2)))))/(length(unique(skipmissing(Ocasião_1)))))*((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1)))))*((sum(D.*J)/(((length(Unidades))-(length(unique(skipmissing(Ocasião_1)))))-1))/(sqrt(Sxm²)*sqrt(Sym²)))^2)))+
        (((((length(Unidades))-(length(unique(skipmissing(Ocasião_1)))))*(((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1)))))/(length(unique(skipmissing(Ocasião_1))))))/
        ((length(unique(skipmissing(Ocasião_2))))-((((length(Unidades))-(length(unique(skipmissing(Ocasião_2)))))/(length(unique(skipmissing(Ocasião_1)))))*
        ((length(Unidades))-(length(unique(skipmissing(Ocasião_1)))))*((sum(D.*J)/(((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1)))))-1))/(sqrt(Sxm²)*sqrt(Sym²)))^2)))*
        (((sum(D.*J)/(((length(Unidades))-(length(unique(skipmissing(Ocasião_1)))))-1))/(sqrt(Sxm²)*sqrt(Sym²)))*(sqrt(Sxm²)/sqrt(Sym²))))
        
        B=((-((length(Unidades))-(length(unique(skipmissing(Ocasião_1)))))*(((length(Unidades))-
        (length(unique(skipmissing(Ocasião_2)))))/(length(unique(skipmissing(Ocasião_1))))))/
        ((length(unique(skipmissing(Ocasião_2))))-((((length(Unidades))-(length(unique(skipmissing(Ocasião_2)))))/
        (length(unique(skipmissing(Ocasião_1)))))*((length(Unidades))-(length(unique(skipmissing(Ocasião_1)))))*((sum(D.*J)/(((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1)))))-1))/(sqrt(Sxm²)*sqrt(Sym²)))^2)))*
        (((sum(D.*J)/(((length(Unidades))-(length(unique(skipmissing(Ocasião_1)))))-1))/(sqrt(Sxm²)*sqrt(Sym²)))*
        (sqrt(Sym²)/sqrt(Sxm²)))-(((length(unique(skipmissing(Ocasião_2))))*(((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1)))))/(length(unique(skipmissing(Ocasião_1))))))/
        ((length(unique(skipmissing(Ocasião_2))))-((((length(Unidades))-(length(unique(skipmissing(Ocasião_2)))))/
        (length(unique(skipmissing(Ocasião_1)))))*((length(Unidades))-(length(unique(skipmissing(Ocasião_1)))))*((sum(D.*J)/(((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1)))))-1))/(sqrt(Sxm²)*sqrt(Sym²)))^2)))

    
        
        #Melhor estimativa da média direta da primeira ocasião
        #Para isso é necessário encontrar os coeficientes "b" "c" dados pelas seguintes equações:
        
        b=((length(unique(skipmissing(Ocasião_2))))*(((length(Unidades))-(length(unique(skipmissing(Ocasião_1)))))/
        (length(unique(skipmissing(Ocasião_1))))))/((length(unique(skipmissing(Ocasião_2))))-
        ((((length(Unidades))-(length(unique(skipmissing(Ocasião_2)))))/(length(unique(skipmissing(Ocasião_1)))))*
        ((length(Unidades))-(length(unique(skipmissing(Ocasião_1)))))*((sum(D.*J)/(((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1)))))-1))/(sqrt(Sxm²)*sqrt(Sym²)))^2))
        
        c=((-((length(Unidades))-(length(unique(skipmissing(Ocasião_1)))))*(((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1)))))/(length(unique(skipmissing(Ocasião_1))))))/
        ((length(unique(skipmissing(Ocasião_2))))-((((length(Unidades))-(length(unique(skipmissing(Ocasião_2)))))/
        (length(unique(skipmissing(Ocasião_1)))))*((length(Unidades))-(length(unique(skipmissing(Ocasião_1)))))*((sum(D.*J)/(((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1)))))-1))/(sqrt(Sxm²)*sqrt(Sym²)))^2)))*(sqrt(Sxm²)/sqrt(Sym²))
        
        #A partir deste ponto, encontra-se a média direta da primeira ocasião
        X=((1-b)*Xu)+(b*Xm)+(c*Ym)-(c*Yn) 
        
        Mudança_crescimento = DataFrame(Variáveis=["Crescimento médio ($(unit)/ha)", "Média direta da primeira ocasião ($(unit)/ha)",
        "Limite inferior do intervalo de confiança para média ($(unit)/ha)", "Limite superior do intervalo de confiança para média ($(unit)/ha)", 
        "Crescimento total estimado ($(unit))", "Limite inferior do intervalo de confiança para o total ($(unit))", 
        "Limite superior do intervalo de confiança para o total ($(unit))", "Erro da amostragem relativo (%)", "Erro da amostragem absoluto ($(unit)/ha)", 
        "Erro padrão ($(unit)/ha)", "Variância da média", "Coeficientes A", "Coeficiente B", "Coeficientes b", "Coeficiente c"], 
        Valores=[(A*Ym)+((1-A)*Yn)+(B*Xm)-((1+B)*Xu), X, ((((a*Xu)-(a*Xm)+(c*Ym)+(1-c)*Yn)-X)-((quantile(TDist((((length(unique(skipmissing(Ocasião_1))))-1)+
        ((length(unique(skipmissing(Ocasião_2))))-1))-1),1-α/2))*(sqrt((A^2)*((sum((skipmissing(Ocasião_2).-(sum(skipmissing(Ocasião_2))/
        (length(unique(skipmissing(Ocasião_2)))))).^2)/(length(unique(skipmissing(Ocasião_2))))-1)/((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1))))))+((1-A)^2)*((sum((skipmissing(Ocasião_2).-(sum(skipmissing(Ocasião_2))/
        (length(unique(skipmissing(Ocasião_2)))))).^2)/(length(unique(skipmissing(Ocasião_2))))-1)/((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1))))))+(B^2)*((sum((skipmissing(Ocasião_1).-(sum(skipmissing(Ocasião_1))/
        (length(unique(skipmissing(Ocasião_1)))))).^2)/(length(unique(skipmissing(Ocasião_1))))-1)/((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1))))))+((1+B)^2)*((sum((skipmissing(Ocasião_1).-(sum(skipmissing(Ocasião_1))/
        (length(unique(skipmissing(Ocasião_1)))))).^2)/(length(unique(skipmissing(Ocasião_1))))-1)/((length(Unidades))-
        (length(unique(skipmissing(Ocasião_2))))))+2*A*B*((sum(D.*J)/(((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1)))))-1))/(sqrt(Sxm²)*sqrt(Sym²)))*((sqrt(Sxm²)*sqrt(Sym²))/((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1)))))))))), (((quantile(TDist((((length(unique(skipmissing(Ocasião_1))))-1)+((length(unique(skipmissing(Ocasião_2))))-1))-1),1-α/2))*(sqrt((A^2)*
        ((sum((skipmissing(Ocasião_2).-(sum(skipmissing(Ocasião_2))/(length(unique(skipmissing(Ocasião_2)))))).^2)/
        (length(unique(skipmissing(Ocasião_2))))-1)/((length(Unidades))-(length(unique(skipmissing(Ocasião_1))))))+((1-A)^2)*
        ((sum((skipmissing(Ocasião_2).-(sum(skipmissing(Ocasião_2))/(length(unique(skipmissing(Ocasião_2)))))).^2)/
        (length(unique(skipmissing(Ocasião_2))))-1)/((length(Unidades))-(length(unique(skipmissing(Ocasião_1))))))+(B^2)*
        ((sum((skipmissing(Ocasião_1).-(sum(skipmissing(Ocasião_1))/(length(unique(skipmissing(Ocasião_1)))))).^2)/
        (length(unique(skipmissing(Ocasião_1))))-1)/((length(Unidades))-(length(unique(skipmissing(Ocasião_1))))))+((1+B)^2)*
        ((sum((skipmissing(Ocasião_1).-(sum(skipmissing(Ocasião_1))/(length(unique(skipmissing(Ocasião_1)))))).^2)/
        (length(unique(skipmissing(Ocasião_1))))-1)/((length(Unidades))-(length(unique(skipmissing(Ocasião_2))))))+2*A*B*((sum(D.*J)/(((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1)))))-1))/(sqrt(Sxm²)*sqrt(Sym²)))*((sqrt(Sxm²)*sqrt(Sym²))/((length(Unidades))-(length(unique(skipmissing(Ocasião_1)))))))))/
        (((a*Xu)-(a*Xm)+(c*Ym)+(1-c)*Yn)-X))*100, N*(((a*Xu)-(a*Xm)+(c*Ym)+(1-c)*Yn)-X), ((N*(((a*Xu)-(a*Xm)+(c*Ym)+(1-c)*Yn)-X))-N*((quantile(TDist((((length(unique(skipmissing(Ocasião_1))))-1)+
        ((length(unique(skipmissing(Ocasião_2))))-1))-1),1-α/2))*(sqrt((A^2)*((sum((skipmissing(Ocasião_2).-(sum(skipmissing(Ocasião_2))/
        (length(unique(skipmissing(Ocasião_2)))))).^2)/(length(unique(skipmissing(Ocasião_2))))-1)/((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1))))))+((1-A)^2)*((sum((skipmissing(Ocasião_2).-(sum(skipmissing(Ocasião_2))/
        (length(unique(skipmissing(Ocasião_2)))))).^2)/(length(unique(skipmissing(Ocasião_2))))-1)/((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1))))))+(B^2)*((sum((skipmissing(Ocasião_1).-(sum(skipmissing(Ocasião_1))/
        (length(unique(skipmissing(Ocasião_1)))))).^2)/(length(unique(skipmissing(Ocasião_1))))-1)/((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1))))))+((1+B)^2)*((sum((skipmissing(Ocasião_1).-(sum(skipmissing(Ocasião_1))/
        (length(unique(skipmissing(Ocasião_1)))))).^2)/(length(unique(skipmissing(Ocasião_1))))-1)/((length(Unidades))-
        (length(unique(skipmissing(Ocasião_2))))))+2*A*B*((sum(D.*J)/(((length(Unidades))-(length(unique(skipmissing(Ocasião_1)))))-1))/(sqrt(Sxm²)*sqrt(Sym²)))*
        ((sqrt(Sxm²)*sqrt(Sym²))/((length(Unidades))-(length(unique(skipmissing(Ocasião_1)))))))))), ((N*(((a*Xu)-(a*Xm)+(c*Ym)+(1-c)*Yn)-X))+N*((quantile(TDist((((length(unique(skipmissing(Ocasião_1))))-1)+
        ((length(unique(skipmissing(Ocasião_2))))-1))-1),1-α/2))*(sqrt((A^2)*((sum((skipmissing(Ocasião_2).-(sum(skipmissing(Ocasião_2))/
        (length(unique(skipmissing(Ocasião_2)))))).^2)/(length(unique(skipmissing(Ocasião_2))))-1)/((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1))))))+((1-A)^2)*((sum((skipmissing(Ocasião_2).-(sum(skipmissing(Ocasião_2))/
        (length(unique(skipmissing(Ocasião_2)))))).^2)/(length(unique(skipmissing(Ocasião_2))))-1)/((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1))))))+(B^2)*((sum((skipmissing(Ocasião_1).-(sum(skipmissing(Ocasião_1))/
        (length(unique(skipmissing(Ocasião_1)))))).^2)/(length(unique(skipmissing(Ocasião_1))))-1)/((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1))))))+((1+B)^2)*((sum((skipmissing(Ocasião_1).-(sum(skipmissing(Ocasião_1))/
        (length(unique(skipmissing(Ocasião_1)))))).^2)/(length(unique(skipmissing(Ocasião_1))))-1)/((length(Unidades))-
        (length(unique(skipmissing(Ocasião_2))))))+2*A*B*((sum(D.*J)/(((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1)))))-1))/(sqrt(Sxm²)*sqrt(Sym²)))*((sqrt(Sxm²)*sqrt(Sym²))/((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1)))))))))), (((quantile(TDist((((length(unique(skipmissing(Ocasião_1))))-1)+((length(unique(skipmissing(Ocasião_2))))-1))-1),1-α/2))*(sqrt((A^2)*
        ((sum((skipmissing(Ocasião_2).-(sum(skipmissing(Ocasião_2))/(length(unique(skipmissing(Ocasião_2)))))).^2)/
        (length(unique(skipmissing(Ocasião_2))))-1)/((length(Unidades))-(length(unique(skipmissing(Ocasião_1))))))+((1-A)^2)*
        ((sum((skipmissing(Ocasião_2).-(sum(skipmissing(Ocasião_2))/(length(unique(skipmissing(Ocasião_2)))))).^2)/
        (length(unique(skipmissing(Ocasião_2))))-1)/((length(Unidades))-(length(unique(skipmissing(Ocasião_1))))))+(B^2)*
        ((sum((skipmissing(Ocasião_1).-(sum(skipmissing(Ocasião_1))/(length(unique(skipmissing(Ocasião_1)))))).^2)/
        (length(unique(skipmissing(Ocasião_1))))-1)/((length(Unidades))-(length(unique(skipmissing(Ocasião_1))))))+((1+B)^2)*
        ((sum((skipmissing(Ocasião_1).-(sum(skipmissing(Ocasião_1))/(length(unique(skipmissing(Ocasião_1)))))).^2)/
        (length(unique(skipmissing(Ocasião_1))))-1)/((length(Unidades))-(length(unique(skipmissing(Ocasião_2))))))+2*A*B*((sum(D.*J)/(((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1)))))-1))/(sqrt(Sxm²)*sqrt(Sym²)))*((sqrt(Sxm²)*sqrt(Sym²))/((length(Unidades))-(length(unique(skipmissing(Ocasião_1)))))))))/
        (((a*Xu)-(a*Xm)+(c*Ym)+(1-c)*Yn)-X))*100, (quantile(TDist((((length(unique(skipmissing(Ocasião_1))))-1)+((length(unique(skipmissing(Ocasião_2))))-1))-1),1-α/2))*(sqrt((A^2)*
        ((sum((skipmissing(Ocasião_2).-(sum(skipmissing(Ocasião_2))/(length(unique(skipmissing(Ocasião_2)))))).^2)/
        (length(unique(skipmissing(Ocasião_2))))-1)/((length(Unidades))-(length(unique(skipmissing(Ocasião_1))))))+((1-A)^2)*
        ((sum((skipmissing(Ocasião_2).-(sum(skipmissing(Ocasião_2))/(length(unique(skipmissing(Ocasião_2)))))).^2)/
        (length(unique(skipmissing(Ocasião_2))))-1)/((length(Unidades))-(length(unique(skipmissing(Ocasião_1))))))+(B^2)*
        ((sum((skipmissing(Ocasião_1).-(sum(skipmissing(Ocasião_1))/(length(unique(skipmissing(Ocasião_1)))))).^2)/
        (length(unique(skipmissing(Ocasião_1))))-1)/((length(Unidades))-(length(unique(skipmissing(Ocasião_1))))))+((1+B)^2)*
        ((sum((skipmissing(Ocasião_1).-(sum(skipmissing(Ocasião_1))/(length(unique(skipmissing(Ocasião_1)))))).^2)/
        (length(unique(skipmissing(Ocasião_1))))-1)/((length(Unidades))-(length(unique(skipmissing(Ocasião_2))))))+2*A*B*((sum(D.*J)/(((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1)))))-1))/(sqrt(Sxm²)*sqrt(Sym²)))*((sqrt(Sxm²)*sqrt(Sym²))/
        ((length(Unidades))-(length(unique(skipmissing(Ocasião_1)))))))), sqrt((A^2)*((sum((skipmissing(Ocasião_2).-(sum(skipmissing(Ocasião_2))/(length(unique(skipmissing(Ocasião_2)))))).^2)/
        (length(unique(skipmissing(Ocasião_2))))-1)/((length(Unidades))-(length(unique(skipmissing(Ocasião_1))))))+((1-A)^2)*
        ((sum((skipmissing(Ocasião_2).-(sum(skipmissing(Ocasião_2))/(length(unique(skipmissing(Ocasião_2)))))).^2)/
        (length(unique(skipmissing(Ocasião_2))))-1)/((length(Unidades))-(length(unique(skipmissing(Ocasião_1))))))+(B^2)*
        ((sum((skipmissing(Ocasião_1).-(sum(skipmissing(Ocasião_1))/(length(unique(skipmissing(Ocasião_1)))))).^2)/
        (length(unique(skipmissing(Ocasião_1))))-1)/((length(Unidades))-(length(unique(skipmissing(Ocasião_1))))))+((1+B)^2)*
        ((sum((skipmissing(Ocasião_1).-(sum(skipmissing(Ocasião_1))/(length(unique(skipmissing(Ocasião_1)))))).^2)/
        (length(unique(skipmissing(Ocasião_1))))-1)/((length(Unidades))-(length(unique(skipmissing(Ocasião_2))))))+2*A*B*((sum(D.*J)/(((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1)))))-1))/(sqrt(Sxm²)*sqrt(Sym²)))*((sqrt(Sxm²)*sqrt(Sym²))/
        ((length(Unidades))-(length(unique(skipmissing(Ocasião_1))))))), ((sum((skipmissing(Ocasião_1).-(sum(skipmissing(Ocasião_1))/(length(unique(skipmissing(Ocasião_1)))))).^2)/
        (length(unique(skipmissing(Ocasião_1))))-1)/(length(unique(skipmissing(Ocasião_1)))))*(1-((((length(Unidades))-(length(unique(skipmissing(Ocasião_1)))))*((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1)))))*((sum(D.*J)/(((length(Unidades))-(length(unique(skipmissing(Ocasião_1)))))-1))/(sqrt(Sxm²)*sqrt(Sym²)))^2)/
        (((length(unique(skipmissing(Ocasião_1))))*(length(unique(skipmissing(Ocasião_2)))))-(((length(Unidades))-(length(unique(skipmissing(Ocasião_2)))))*((length(Unidades))-
        (length(unique(skipmissing(Ocasião_1)))))*((sum(D.*J)/(((length(Unidades))-(length(unique(skipmissing(Ocasião_1)))))-1))/
        (sqrt(Sxm²)*sqrt(Sym²)))^2)))), A, B, b, c]) #Tabela de resultados
        
        Resultados = [ARP, Informações_do_inventário, Primeira_ocasião, Segunda_ocasião, Mudança_crescimento]
        
        return [Resultados]
    end

    function ART(Dados, AreaParc, N1, N2, n_unidade, Ocasiao_1, Ocasiao_2, α, unit)

        AreaParc = Float64(Meta.parse(AreaParc))
        α = Float64(Meta.parse(α))
        N1 = Int64(Meta.parse(N1))
        N2 = Int64(Meta.parse(N2))

        Conversor=1/AreaParc

        ###Primeira ocasião####
        Unidades = Dados[!,n_unidade]
        Ocasião_1 = (Conversor.*Dados[!, Ocasiao_1])
        Ocasião_2 = (Conversor.*Dados[!, Ocasiao_2])
        Independente = DataFrame(Unidades = Unidades, Ocasião_1 = Ocasião_1, Ocasião_2 = Ocasião_2)
    
        Primeira_ocasião = DataFrame(Variáveis=["Média ($(unit)/ha)", "Limite inferior do intervalo de confiança para média ($(unit)/ha)", 
        "Limite superior do intervalo de confiança para média ($(unit)/ha)", "Total da população ($(unit))", "Limite inferior do intervalo de confiança para o total ($(unit))", 
        "Limite superior do intervalo de confiança para o total ($(unit))", "Área da população (ha)", "Erro da amostragem relativo (%)", 
        "Erro da amostragem absoluto ($(unit)/ha)", "Erro padrão ($(unit)/ha)", "Desvio padrão ($(unit)/ha)", "Variância ($(unit)/ha)²", "Variância da média ($(unit)/ha)²", 
        "Número total de unidades", "Nível de de significância (α)"], Valores=[mean(Ocasião_1), ((mean(Ocasião_1))-((quantile(TDist((length(Ocasião_1))-1),1-α/2))*
        (sqrt(((sum((Ocasião_1.-(mean(Ocasião_1))).^2)/(length(Ocasião_1))-1)/(length(Ocasião_1)))*(1-((length(Ocasião_1))/N1)))))),
        ((mean(Ocasião_1))+((quantile(TDist((length(Ocasião_1))-1),1-α/2))*(sqrt(((sum((Ocasião_1.-(mean(Ocasião_1))).^2)/
        (length(Ocasião_1))-1)/(length(Ocasião_1)))*(1-((length(Ocasião_1))/N1)))))), (N1*(mean(Ocasião_1))/Conversor), (((N1*(mean(Ocasião_1)))-
        N1*((quantile(TDist((length(Ocasião_1))-1),1-α/2))*(sqrt(((sum((Ocasião_1.-(mean(Ocasião_1))).^2)/(length(Ocasião_1))-1)/
        (length(Ocasião_1)))*(1-((length(Ocasião_1))/N1))))))/Conversor), (((N1*(mean(Ocasião_1)))+N1*((quantile(TDist((length(Ocasião_1))-1),1-α/2))*
        (sqrt(((sum((Ocasião_1.-(mean(Ocasião_1))).^2)/(length(Ocasião_1))-1)/(length(Ocasião_1)))*(1-((length(Ocasião_1))/N1))))))/Conversor), N1, 
        (((quantile(TDist((length(Ocasião_1))-1),1-α/2))*(sqrt(((sum((Ocasião_1.-(mean(Ocasião_1))).^2)/(length(Ocasião_1))-1)/(length(Ocasião_1)))*
        (1-((length(Ocasião_1))/N1)))))/(mean(Ocasião_1)))*100, (quantile(TDist((length(Ocasião_1))-1),1-α/2))*(sqrt(((sum((Ocasião_1.-(mean(Ocasião_1))).^2)/
        (length(Ocasião_1))-1)/(length(Ocasião_1)))*(1-((length(Ocasião_1))/N1)))), sqrt(((sum((Ocasião_1.-(mean(Ocasião_1))).^2)/(length(Ocasião_1))-1)/(length(Ocasião_1)))*
        (1-((length(Ocasião_1))/N1))), sqrt(sum((Ocasião_1.-(mean(Ocasião_1))).^2)/(length(Ocasião_1))-1), sum((Ocasião_1.-(mean(Ocasião_1))).^2)/(length(Ocasião_1))-1, 
        ((sum((Ocasião_1.-(mean(Ocasião_1))).^2)/(length(Ocasião_1))-1)/(length(Ocasião_1)))*(1-((length(Ocasião_1))/N1)), length(Ocasião_1), α]) #Tabela de resultados
        
        Segunda_ocasião = DataFrame(Variáveis=["Média ($(unit)/ha)", "Limite inferior do intervalo de confiança para média ($(unit)/ha)", 
        "Limite superior do intervalo de confiança para média ($(unit)/ha)", "Total da população ($(unit))", "Limite inferior do intervalo de confiança para o total ($(unit))", 
        "Limite superior do intervalo de confiança para o total ($(unit))", "Área da população (ha)", "Erro da amostragem relativo (%)", 
        "Erro da amostragem absoluto ($(unit)/ha)", "Erro padrão ($(unit)/ha)", "Desvio padrão ($(unit)/ha)", "Variância ($(unit)/ha)²", "Variância da média ($(unit)/ha)²", 
        "Número total de unidades", "Nível de de significância (α)"], Valores=[mean(Ocasião_2), ((mean(Ocasião_2))-((quantile(TDist((length(Ocasião_2))-1),1-α/2))*
        (sqrt(((sum((Ocasião_2.-(mean(Ocasião_2))).^2)/(length(Ocasião_2))-1)/(length(Ocasião_2)))*(1-((length(Ocasião_2))/N2)))))), 
        ((mean(Ocasião_2))+((quantile(TDist((length(Ocasião_2))-1),1-α/2))*(sqrt(((sum((Ocasião_2.-(mean(Ocasião_2))).^2)/
        (length(Ocasião_2))-1)/(length(Ocasião_2)))*(1-((length(Ocasião_2))/N2)))))), (N2*(mean(Ocasião_2))/Conversor), 
        (((N2*(mean(Ocasião_2)))-N2*((quantile(TDist((length(Ocasião_2))-1),1-α/2))*(sqrt(((sum((Ocasião_2.-(mean(Ocasião_2))).^2)/
        (length(Ocasião_2))-1)/(length(Ocasião_2)))*(1-((length(Ocasião_2))/N2))))))/Conversor), ((N2*((mean(Ocasião_2)))+N2*
        ((quantile(TDist((length(Ocasião_2))-1),1-α/2))*(sqrt(((sum((Ocasião_2.-(mean(Ocasião_2))).^2)/(length(Ocasião_2))-1)/(length(Ocasião_2)))*
        (1-((length(Ocasião_2))/N2))))))/Conversor), N2, (((quantile(TDist((length(Ocasião_2))-1),1-α/2))*(sqrt(((sum((Ocasião_2.-(mean(Ocasião_2))).^2)/
        (length(Ocasião_2))-1)/(length(Ocasião_2)))*(1-((length(Ocasião_2))/N2)))))/(mean(Ocasião_2)))*100, 
        (quantile(TDist((length(Ocasião_2))-1),1-α/2))*(sqrt(((sum((Ocasião_2.-(mean(Ocasião_2))).^2)/
        (length(Ocasião_2))-1)/(length(Ocasião_2)))*(1-((length(Ocasião_2))/N2)))), sqrt(((sum((Ocasião_2.-(mean(Ocasião_2))).^2)/(length(Ocasião_2))-1)/
        (length(Ocasião_2)))*(1-((length(Ocasião_2))/N2))), sqrt(sum((Ocasião_2.-(mean(Ocasião_2))).^2)/(length(Ocasião_2))-1), 
        sum((Ocasião_2.-(mean(Ocasião_2))).^2)/(length(Ocasião_2))-1, ((sum((Ocasião_2.-(mean(Ocasião_2))).^2)/(length(Ocasião_2))-1)/
        (length(Ocasião_2)))*(1-((length(Ocasião_2))/N2)), length(Ocasião_2), α]) #Tabela de resultados
        
        Mudança_crescimento = DataFrame(Variáveis=["Crescimento médio ($(unit)/ha)", "Limite inferior do intervalo de confiança para média ($(unit)/ha)", 
        "Limite superior do intervalo de confiança para média ($(unit)/ha)", "Crescimento total estimado ($(unit))", "Limite inferior do intervalo de confiança para o total ($(unit))", 
        "Limite superior do intervalo de confiança para o total ($(unit))", "Área da população (ha)", "Erro da amostragem relativo (%)", 
        "Erro da amostragem absoluto ($(unit)/ha)", "Erro padrão ($(unit)/ha)", "Variância da média ($(unit)/ha)²"], Valores=[mean(Ocasião_2)-mean(Ocasião_1), 
        ((mean(Ocasião_2)-mean(Ocasião_1))-quantile(TDist((((length(Ocasião_1))-1)+((length(Ocasião_2))-1))-1),1-α/2)*
        (sqrt(((((sum((Ocasião_1.-(mean(Ocasião_1))).^2)/(length(Ocasião_1))-1)/(length(Ocasião_1)))*
        (1-((length(Ocasião_1))/N1)))+(((sum((Ocasião_2.-(mean(Ocasião_2))).^2)/(length(Ocasião_2))-1)/(length(Ocasião_2)))*(1-((length(Ocasião_2))/N2))))-
        (2*((sum((Ocasião_1.-(mean(Ocasião_1))).*(Ocasião_2.-(mean(Ocasião_2))))/(length(Dados.n)-1))/
        (length(Dados.n))))))), ((mean(Ocasião_2)-mean(Ocasião_1))+quantile(TDist((((length(Ocasião_1))-1)+((length(Ocasião_2))-1))-1),1-α/2)*
        (sqrt(((((sum((Ocasião_1.-(mean(Ocasião_1))).^2)/(length(Ocasião_1))-1)/(length(Ocasião_1)))*
        (1-((length(Ocasião_1))/N1)))+(((sum((Ocasião_2.-(mean(Ocasião_2))).^2)/(length(Ocasião_2))-1)/(length(Ocasião_2)))*(1-((length(Ocasião_2))/N2))))-
        (2*((sum((Ocasião_1.-(mean(Ocasião_1))).*(Ocasião_2.-(mean(Ocasião_2))))/(length(Dados.n)-1))/
        (length(Dados.n))))))), N2*(mean(Ocasião_2)-mean(Ocasião_1)), ((N2*(mean(Ocasião_2)-mean(Ocasião_1)))-N2*quantile(TDist((((length(Ocasião_1))-1)+
        ((length(Ocasião_2))-1))-1),1-α/2)*(sqrt(((((sum((Ocasião_1.-(mean(Ocasião_1))).^2)/
        (length(Ocasião_1))-1)/(length(Ocasião_1)))*(1-((length(Ocasião_1))/N1)))+(((sum((Ocasião_2.-(mean(Ocasião_2))).^2)/(length(Ocasião_2))-1)/(length(Ocasião_2)))*
        (1-((length(Ocasião_2))/N2))))-(2*((sum((Ocasião_1.-(mean(Ocasião_1))).*(Ocasião_2.-(mean(Ocasião_2))))/(length(Dados.n)-1))/
        (length(Dados.n))))))), ((N2*(mean(Ocasião_2)-mean(Ocasião_1)))+N2*quantile(TDist((((length(Ocasião_1))-1)+
        ((length(Ocasião_2))-1))-1),1-α/2)*(sqrt(((((sum((Ocasião_1.-(mean(Ocasião_1))).^2)/
        (length(Ocasião_1))-1)/(length(Ocasião_1)))*(1-((length(Ocasião_1))/N1)))+(((sum((Ocasião_2.-(mean(Ocasião_2))).^2)/(length(Ocasião_2))-1)/(length(Ocasião_2)))*
        (1-((length(Ocasião_2))/N2))))-(2*((sum((Ocasião_1.-(mean(Ocasião_1))).*(Ocasião_2.-(mean(Ocasião_2))))/(length(Dados.n)-1))/
        (length(Dados.n))))))), N2, (quantile(TDist((((length(Ocasião_1))-1)+((length(Ocasião_2))-1))-1),1-α/2)*
        (sqrt(((((sum((Ocasião_1.-(mean(Ocasião_1))).^2)/(length(Ocasião_1))-1)/(length(Ocasião_1)))*
        (1-((length(Ocasião_1))/N1)))+(((sum((Ocasião_2.-(mean(Ocasião_2))).^2)/(length(Ocasião_2))-1)/
        (length(Ocasião_2)))*(1-((length(Ocasião_2))/N2))))-(2*((sum((Ocasião_1.-(mean(Ocasião_1))).*(Ocasião_2.-(mean(Ocasião_2))))/(length(Dados.n)-1))/
        (length(Dados.n))))))/(mean(Ocasião_2)-mean(Ocasião_1)))*100, quantile(TDist((((length(Ocasião_1))-1)+((length(Ocasião_2))-1))-1),1-α/2)*
        (sqrt(((((sum((Ocasião_1.-(mean(Ocasião_1))).^2)/(length(Ocasião_1))-1)/(length(Ocasião_1)))*
        (1-((length(Ocasião_1))/N1)))+(((sum((Ocasião_2.-(mean(Ocasião_2))).^2)/(length(Ocasião_2))-1)/
        (length(Ocasião_2)))*(1-((length(Ocasião_2))/N2))))-(2*((sum((Ocasião_1.-(mean(Ocasião_1))).*(Ocasião_2.-(mean(Ocasião_2))))/
        (length(Dados.n)-1))/(length(Dados.n)))))), sqrt(((((sum((Ocasião_1.-(mean(Ocasião_1))).^2)/(length(Ocasião_1))-1)/(length(Ocasião_1)))*
        (1-((length(Ocasião_1))/N1)))+(((sum((Ocasião_2.-(mean(Ocasião_2))).^2)/(length(Ocasião_2))-1)/
        (length(Ocasião_2)))*(1-((length(Ocasião_2))/N2))))-(2*((sum((Ocasião_1.-(mean(Ocasião_1))).*(Ocasião_2.-(mean(Ocasião_2))))/(length(Dados.n)-1))/
        (length(Dados.n))))), ((((sum((Ocasião_1.-(mean(Ocasião_1))).^2)/(length(Ocasião_1))-1)/(length(Ocasião_1)))*
        (1-((length(Ocasião_1))/N1)))+(((sum((Ocasião_2.-(mean(Ocasião_2))).^2)/(length(Ocasião_2))-1)/
        (length(Ocasião_2)))*(1-((length(Ocasião_2))/N2))))-(2*((sum((Ocasião_1.-(mean(Ocasião_1))).*(Ocasião_2.-(mean(Ocasião_2))))/(length(Dados.n)-1))/
        (length(Dados.n))))]) #Tabela de resultados   

        Resultados = [Independente, Primeira_ocasião, Segunda_ocasião, Mudança_crescimento]
        
        return Mudança_crescimento
    end

    function CONGL(Dados, AreaParc, Area, EAR, α, unit) #Determinar função

        Area = Float64(Meta.parse(Area))
        AreaParc = Float64(Meta.parse(AreaParc))
        α = Float64(Meta.parse(α))
        EAR = Float64(Meta.parse(EAR))

        N = (Area*10000)

        Conversor = 10000/AreaParc

        Conjunto_de_dados = (Conversor.*Dados)

        #Tabela com estatítica descritiva por unidades/blocos secundários
        Tabela=transform(Conjunto_de_dados, AsTable(:) .=> ByRow.([I -> count(!ismissing, I), sum, mean, var]).=>[:n, :Soma, :Média, :Variância])
        
        if (quantile(TDist(length(Tabela.n)-1),1-α/2)*sqrt((((N-length(Tabela.n))/N))*(((sum(first(unique(Tabela.n)).*
            (Tabela.Média.-sum(Tabela.Média)/(length(Tabela.n))).^2)/(length(Tabela.n)-1).-sum(Tabela.Variância.*(first(unique(Tabela.n)).-1))/
            (length(Tabela.n)*(first(unique(Tabela.n)).-1)))./first(unique(Tabela.n))/length(Tabela.n))+
            (sum(Tabela.Variância.*(first(unique(Tabela.n)).-1))/(length(Tabela.n)*(first(unique(Tabela.n)).-1))./
            (length(Tabela.n)*first(unique(Tabela.n))))))/(sum(Tabela.Média)/(length(Tabela.n))))*100 > EAR
            Observação = "Diante do exposto, conclui-se que os resultados obtidos na amostragem não satisfazem as exigências de precisão estabelecidas para o inventário, ou seja, um erro de amostragem máximo de $(Int(EAR))% da média  para confiabilidade designada. \n\nO erro estimado foi maior que o limite fixado, sendo recomendado incluir mais unidades amostrais no inventário."
        elseif (quantile(TDist(length(Tabela.n)-1),1-α/2)*sqrt((((N-length(Tabela.n))/N))*(((sum(first(unique(Tabela.n)).*
            (Tabela.Média.-sum(Tabela.Média)/(length(Tabela.n))).^2)/(length(Tabela.n)-1).-sum(Tabela.Variância.*(first(unique(Tabela.n)).-1))/
            (length(Tabela.n)*(first(unique(Tabela.n)).-1)))./first(unique(Tabela.n))/length(Tabela.n))+
            (sum(Tabela.Variância.*(first(unique(Tabela.n)).-1))/(length(Tabela.n)*(first(unique(Tabela.n)).-1))./
            (length(Tabela.n)*first(unique(Tabela.n))))))/(sum(Tabela.Média)/(length(Tabela.n))))*100 ≤ EAR
            Observação  = "Diante do exposto, conclui-se que os resultados obtidos na amostragem satisfazem as exigências de precisão estabelecidas para o inventário, ou seja, um erro de amostragem máximo de ±$(Int(EAR))% da média para confiabilidade designada. \n\nO erro estimado foi menor que o limite fixado, assim as unidades amostrais são suficientes para o inventário."
        end   

        Resultados = DataFrame(Variáveis=["Média ($(unit)/ha)", "Limite inferior do intervalo de confiança para média ($(unit)/ha)", 
        "Limite superior do intervalo de confiança para média ($(unit)/ha)", "Total da população ($(unit))", "Limite inferior do intervalo de confiança para o total ($(unit))", 
        "Limite superior do intervalo de confiança para o total ($(unit))", "Área da população (ha)", "Erro padrão relativo (%)", "Erro da amostragem absoluto ($(unit)/ha)", "Erro padrão ($(unit)/ha)",
        "Variância dentro dos conglomerados ($(unit)/ha)²", "Variância entre conglomerados ($(unit)/ha)²", "Variância da população por subunidade ($(unit)/ha)²", 
        "Variância da população total ($(unit)/ha)²", "Variância da média ($(unit)/ha)²", "Coeficiente de correlação intraconglomerados", "Tamanho da amostra",
        "Limite do erro de amostragem requerido", "Número de unidades primárias", "Número de unidades secundarias", "Nível de significância (α)", "Observação"], 
        Valores=[sum(Tabela.Média)/(length(Tabela.n)), ((sum(Tabela.Média)/(length(Tabela.n)))-(quantile(TDist(length(Tabela.n)-1),1-α/2)*sqrt((((N-length(Tabela.n))/N))*
        (((sum(first(unique(Tabela.n)).*(Tabela.Média.-sum(Tabela.Média)/(length(Tabela.n))).^2)/(length(Tabela.n)-1).-sum(Tabela.Variância.*(first(unique(Tabela.n)).-1))/(length(Tabela.n)*
        (first(unique(Tabela.n)).-1)))./first(unique(Tabela.n))/length(Tabela.n))+(sum(Tabela.Variância.*(first(unique(Tabela.n)).-1))/(length(Tabela.n)*
        (first(unique(Tabela.n)).-1))./(length(Tabela.n)*first(unique(Tabela.n))))))/(sum(Tabela.Média)/(length(Tabela.n))))*100), (sum(Tabela.Média)/(length(Tabela.n)))+(quantile(TDist(length(Tabela.n)-1),1-α/2)*sqrt((((N-length(Tabela.n))/N))*
        (((sum(first(unique(Tabela.n)).*(Tabela.Média.-sum(Tabela.Média)/(length(Tabela.n))).^2)/(length(Tabela.n)-1).-sum(Tabela.Variância.*(first(unique(Tabela.n)).-1))/(length(Tabela.n)*
        (first(unique(Tabela.n)).-1)))./first(unique(Tabela.n))/length(Tabela.n))+(sum(Tabela.Variância.*(first(unique(Tabela.n)).-1))/
        (length(Tabela.n)*(first(unique(Tabela.n)).-1))./(length(Tabela.n)*first(unique(Tabela.n))))))), ((N*(first(unique(Tabela.n)))*
        (sum(Tabela.Média)/(length(Tabela.n))))/Conversor), (((N*(first(unique(Tabela.n)))*(sum(Tabela.Média)/(length(Tabela.n))))-
        ((N*first(unique(Tabela.n)))*quantile(TDist(length(Tabela.n)-1),1-α/2)*sqrt((((N-length(Tabela.n))/N))*
        (((sum(first(unique(Tabela.n)).*(Tabela.Média.-sum(Tabela.Média)/(length(Tabela.n))).^2)/
        (length(Tabela.n)-1).-sum(Tabela.Variância.*(first(unique(Tabela.n)).-1))/(length(Tabela.n)*
        (first(unique(Tabela.n)).-1)))./first(unique(Tabela.n))/length(Tabela.n))+(sum(Tabela.Variância.*(first(unique(Tabela.n)).-1))/
        (length(Tabela.n)*(first(unique(Tabela.n)).-1))./(length(Tabela.n)*first(unique(Tabela.n))))))))/Conversor), (((N*(first(unique(Tabela.n)))*
        (sum(Tabela.Média)/(length(Tabela.n))))+((N*first(unique(Tabela.n)))*quantile(TDist(length(Tabela.n)-1),1-α/2)*sqrt((((N-length(Tabela.n))/N))*
        (((sum(first(unique(Tabela.n)).*(Tabela.Média.-sum(Tabela.Média)/(length(Tabela.n))).^2)/
        (length(Tabela.n)-1).-sum(Tabela.Variância.*(first(unique(Tabela.n)).-1))/(length(Tabela.n)*
        (first(unique(Tabela.n)).-1)))./first(unique(Tabela.n))/length(Tabela.n))+(sum(Tabela.Variância.*(first(unique(Tabela.n)).-1))/
        (length(Tabela.n)*(first(unique(Tabela.n)).-1))./(length(Tabela.n)*first(unique(Tabela.n))))))))/Conversor), Area, 
        (quantile(TDist(length(Tabela.n)-1),1-α/2)*sqrt((((N-length(Tabela.n))/N))*(((sum(first(unique(Tabela.n)).*
        (Tabela.Média.-sum(Tabela.Média)/(length(Tabela.n))).^2)/(length(Tabela.n)-1).-sum(Tabela.Variância.*(first(unique(Tabela.n)).-1))/
        (length(Tabela.n)*(first(unique(Tabela.n)).-1)))./first(unique(Tabela.n))/length(Tabela.n))+
        (sum(Tabela.Variância.*(first(unique(Tabela.n)).-1))/(length(Tabela.n)*(first(unique(Tabela.n)).-1))./
        (length(Tabela.n)*first(unique(Tabela.n))))))/(sum(Tabela.Média)/(length(Tabela.n))))*100, quantile(TDist(length(Tabela.n)-1),1-α/2)*sqrt((((N-length(Tabela.n))/N))*(((sum(first(unique(Tabela.n)).*
        (Tabela.Média.-sum(Tabela.Média)/(length(Tabela.n))).^2)/(length(Tabela.n)-1).-sum(Tabela.Variância.*(first(unique(Tabela.n)).-1))/
        (length(Tabela.n)*(first(unique(Tabela.n)).-1)))./first(unique(Tabela.n))/length(Tabela.n))+
        (sum(Tabela.Variância.*(first(unique(Tabela.n)).-1))/(length(Tabela.n)*(first(unique(Tabela.n)).-1))./(length(Tabela.n)*first(unique(Tabela.n)))))), 
        sqrt((((N-length(Tabela.n))/N))*(((sum(first(unique(Tabela.n)).*(Tabela.Média.-sum(Tabela.Média)/(length(Tabela.n))).^2)/
        (length(Tabela.n)-1).-sum(Tabela.Variância.*(first(unique(Tabela.n)).-1))/(length(Tabela.n)*(first(unique(Tabela.n)).-1)))./first(unique(Tabela.n))/length(Tabela.n))+
        (sum(Tabela.Variância.*(first(unique(Tabela.n)).-1))/(length(Tabela.n)*(first(unique(Tabela.n)).-1))./(length(Tabela.n)*first(unique(Tabela.n)))))), 
        sum(Tabela.Variância/length(Tabela.n))/Conversor, (((sum(first(unique(Tabela.n)).*(Tabela.Média.-sum(Tabela.Média)/(length(Tabela.n))).^2)/
        (length(Tabela.n)-1).-sum(Tabela.Variância.*(first(unique(Tabela.n)).-1))/(length(Tabela.n)*(first(unique(Tabela.n)).-1)))./first(unique(Tabela.n))))/Conversor, 
        (sum(Tabela.Variância.*(first(unique(Tabela.n)).-1))/(length(Tabela.n)*(first(unique(Tabela.n)).-1))+(sum(first(unique(Tabela.n)).*(Tabela.Média.-sum(Tabela.Média)/
        (length(Tabela.n))).^2)/(length(Tabela.n)-1).-sum(Tabela.Variância.*(first(unique(Tabela.n)).-1))/(length(Tabela.n)*
        (first(unique(Tabela.n)).-1)))./first(unique(Tabela.n)))/Conversor, ((((sum(first(unique(Tabela.n)).*(Tabela.Média.-sum(Tabela.Média)/(length(Tabela.n))).^2)/
        (length(Tabela.n)-1).-sum(Tabela.Variância.*(first(unique(Tabela.n)).-1))/(length(Tabela.n)*
        (first(unique(Tabela.n)).-1)))./first(unique(Tabela.n))))/Conversor)+(sum(Tabela.Variância/length(Tabela.n))/Conversor), (((N-length(Tabela.n))/N))*((((sum(first(unique(Tabela.n)).*(Tabela.Média.-sum(Tabela.Média)/(length(Tabela.n))).^2)/
        (length(Tabela.n)-1).-sum(Tabela.Variância.*(first(unique(Tabela.n)).-1))/(length(Tabela.n)*(first(unique(Tabela.n)).-1)))./first(unique(Tabela.n))))/Conversor)/
        (length(Tabela.n))+(sum(Tabela.Variância.*(first(unique(Tabela.n)).-1))/(length(Tabela.n)*(first(unique(Tabela.n)).-1))/Conversor)./
        (length(Tabela.n)*first(unique(Tabela.n))), (sum(first(unique(Tabela.n)).*(Tabela.Média.-sum(Tabela.Média)/(length(Tabela.n))).^2)/
        (length(Tabela.n)-1).-sum(Tabela.Variância.*(first(unique(Tabela.n)).-1))/(length(Tabela.n)*(first(unique(Tabela.n)).-1)))./first(unique(Tabela.n))/
        ((sum(first(unique(Tabela.n)).*(Tabela.Média.-sum(Tabela.Média)/(length(Tabela.n))).^2)/(length(Tabela.n)-1).-sum(Tabela.Variância.*(first(unique(Tabela.n)).-1))/
        (length(Tabela.n)*(first(unique(Tabela.n)).-1)))./first(unique(Tabela.n))+sum(Tabela.Variância.*(first(unique(Tabela.n)).-1))/
        (length(Tabela.n)*(first(unique(Tabela.n)).-1))), round(((((quantile(TDist(length(Tabela.n)-1),1-α/2))^2)*(sum(Tabela.Variância.*(first(unique(Tabela.n)).-1))/(length(Tabela.n)*
        (first(unique(Tabela.n)).-1))+(sum(first(unique(Tabela.n)).*(Tabela.Média.-sum(Tabela.Média)/(length(Tabela.n))).^2)/(length(Tabela.n)-1).-sum(Tabela.Variância.*(first(unique(Tabela.n)).-1))/
        (length(Tabela.n)*(first(unique(Tabela.n)).-1)))./first(unique(Tabela.n)))))./((((0.1*sum(Tabela.Média)/(length(Tabela.n)))).^2).*(first(unique(Tabela.n))))).*(1+(sum(first(unique(Tabela.n)).*(Tabela.Média.-sum(Tabela.Média)/
        (length(Tabela.n))).^2)/(length(Tabela.n)-1).-sum(Tabela.Variância.*(first(unique(Tabela.n)).-1))/(length(Tabela.n)*(first(unique(Tabela.n)).-1)))./first(unique(Tabela.n))/
        ((sum(first(unique(Tabela.n)).*(Tabela.Média.-sum(Tabela.Média)/(length(Tabela.n))).^2)/(length(Tabela.n)-1).-sum(Tabela.Variância.*(first(unique(Tabela.n)).-1))/
        (length(Tabela.n)*(first(unique(Tabela.n)).-1)))./first(unique(Tabela.n))+sum(Tabela.Variância.*(first(unique(Tabela.n)).-1))/(length(Tabela.n)*(first(unique(Tabela.n)).-1)))'*
        (first(unique(Tabela.n)).-1)), (0.1*sum(Tabela.Média)/(length(Tabela.n))), length(Tabela.n), first(unique(Tabela.n)), α, Observação]) #Tabela de resultados  
        
        #Resultados = [Conjunto_de_dados, Tabela, Resultados]

        return Resultados

    end

    function DE(Dados, AreaParc, Area, EAR, α, M, unit) #Determina a função

        Area = Float64(Meta.parse(Area))
        AreaParc = Float64(Meta.parse(AreaParc))
        α = Float64(Meta.parse(α))
        EAR = Float64(Meta.parse(EAR))
        M = Float64(Meta.parse(M))

        N = Area

        Conversor = 10000/AreaParc

        Conjunto_de_dados = (Conversor.*Dados)

        #Tabela com estatítica descritiva por unidade secundária/bloco
        Tabela=transform(Conjunto_de_dados, AsTable(:) .=> ByRow.([I -> count(!ismissing, I), sum, mean, var]).=>[:n, :Soma, :Média, :Variância])

        if 1-(length(Tabela.n)/N) ≥ 0.98 #f maior ou igual a 0,98 população infinita
            População = "A população avaliada é considerada infinita"   
        elseif 1-(length(Tabela.n)/N) < 0.98 #f menor que 0,98 população finita
            População = "A população avaliada é considerada finita"    
        end
        
        Tamanho_da_amostra = if 1-(length(Tabela.n)/N) ≥ 0.98
             #População infinita. O tamanho da amostra é calculado pela seguinte equação:
             Infinita=((quantile(TDist(length(Tabela.n)-1),1-α/2))^2)*((sum(first(unique(Tabela.n))*
            (Tabela.Média.-(sum(Tabela.Média)/length(Tabela.n))).^2)/
            (length(Tabela.n)-1)-sum(Tabela.Variância*(first(unique(Tabela.n))-1))/(length(Tabela.n)*
            (first(unique(Tabela.n))-1)))/first(unique(Tabela.n))+(sum(Tabela.Variância*(first(unique(Tabela.n))-1))/
            (length(Tabela.n)*(first(unique(Tabela.n))-1))/first(unique(Tabela.n))))/
            (((0.1*(sum(Tabela.Média)/length(Tabela.n))))^2)
            round(Infinita)
        elseif 1-(length(Tabela.n)/N) < 0.98
             #População finita. O tamanho da amaostra é calculado pela seguinte equação:
             Finita=((quantile(TDist(length(Tabela.n)-1),1-α/2))^2)*((sum(first(unique(Tabela.n))*
            (Tabela.Média.-(sum(Tabela.Média)/length(Tabela.n))).^2)/
            (length(Tabela.n)-1)-sum(Tabela.Variância*(first(unique(Tabela.n))-1))/(length(Tabela.n)*
            (first(unique(Tabela.n))-1)))/first(unique(Tabela.n))+(sum(Tabela.Variância*(first(unique(Tabela.n))-1))/
            (length(Tabela.n)*(first(unique(Tabela.n))-1))./first(unique(Tabela.n))))/
            ((((0.1*(sum(Tabela.Média)/length(Tabela.n))))^2)+(1/N)*((quantile(TDist(length(Tabela.n)-1),1-α/2))^2)*
            ((sum(first(unique(Tabela.n))*(Tabela.Média.-(sum(Tabela.Média)/length(Tabela.n))).^2)/
            ((length(Tabela.n))-1)-sum(Tabela.Variância*(first(unique(Tabela.n))-1))/(length(Tabela.n)*
            (first(unique(Tabela.n))-1)))/first(unique(Tabela.n))+(sum(Tabela.Variância*(first(unique(Tabela.n))-1))/
            (length(Tabela.n)*(first(unique(Tabela.n))-1))/M)))
            round(Finita)
        end
        
        if (quantile(TDist(length(Tabela.n)-1),1-α/2)*sqrt((((N-length(Tabela.n))/N)*((sum(first(unique(Tabela.n))*
            (Tabela.Média.-(sum(Tabela.Média)/length(Tabela.n))).^2)/
            (length(Tabela.n)-1)-sum(Tabela.Variância*(first(unique(Tabela.n))-1))/(length(Tabela.n)*
            (first(unique(Tabela.n))-1)))/first(unique(Tabela.n))/length(Tabela.n))+
            ((M-first(unique(Tabela.n)))/M)*(sum(Tabela.Variância*(first(unique(Tabela.n))-1))/(length(Tabela.n)*
            (first(unique(Tabela.n))-1))/(length(Tabela.n)*first(unique(Tabela.n))))))/
            (sum(Tabela.Média)/length(Tabela.n)))*100 > EAR
            Observação = "Diante do exposto, conclui-se que os resultados obtidos na amostragem não satisfazem as exigências de precisão estabelecidas para o inventário, ou seja, um erro de amostragem máximo de ±$(Int(EAR))% da média  para confiabilidade designada. \n\nO erro estimado foi maior que o limite fixado, sendo recomendado incluir mais unidades amostrais no inventário."
        elseif (quantile(TDist(length(Tabela.n)-1),1-α/2)*sqrt((((N-length(Tabela.n))/N)*((sum(first(unique(Tabela.n))*
            (Tabela.Média.-(sum(Tabela.Média)/length(Tabela.n))).^2)/
            (length(Tabela.n)-1)-sum(Tabela.Variância*(first(unique(Tabela.n))-1))/(length(Tabela.n)*
            (first(unique(Tabela.n))-1)))/first(unique(Tabela.n))/length(Tabela.n))+
            ((M-first(unique(Tabela.n)))/M)*(sum(Tabela.Variância*(first(unique(Tabela.n))-1))/(length(Tabela.n)*
            (first(unique(Tabela.n))-1))/(length(Tabela.n)*first(unique(Tabela.n))))))/
            (sum(Tabela.Média)/length(Tabela.n)))*100 ≤ EAR
            Observação  = "Diante do exposto, conclui-se que os resultados obtidos na amostragem satisfazem as exigências de precisão estabelecidas para o inventário, ou seja, um erro de amostragem máximo de ±$(Int(EAR))% da média para confiabilidade designada. \n\nO erro estimado foi menor que o limite fixado, assim as unidades amostrais são suficientes para o inventário."
        end  

        Resultados = DataFrame(Variânciaiáveis=["Média ($(unit)/ha)", "Limite inferior do intervalo de confiança para média ($(unit)/ha)", 
        "Limite superior do intervalo de confiança para média ($(unit)/ha)", "Total da população ($(unit))", "Limite inferior do intervalo de confiança para o total ($(unit))", 
        "Limite superior do intervalo de confiança para o total ($(unit))", "Área da população (ha)", "Erro da amostragem relativo (%)", 
        "Erro da amostragem absoluto ($(unit)/ha)", "Erro padrão ($(unit)/ha)", "Variância dentro das unidades ($(unit)/ha)²", "Variância entre unidades ($(unit)/ha)²", 
        "Estimativa da Variância ($(unit)/ha)²", "Variância da média da população ($(unit)/ha)²", "Limite do erro de amostragem requerido", 
        "Fator de correção", "Tamanho da amostra", "Número total de unidades secundárias por unidade primária", "Número potencial de unidades primárias", 
        "Número de unidades primárias", "Número de unidades secundárias", "Nível de significância (α)"], Valores = round.(([(sum(Tabela.Média)/length(Tabela.n)), 
        ((sum(Tabela.Média)/length(Tabela.n))-quantile(TDist(length(Tabela.n)-1),1-α/2)*sqrt((((N-length(Tabela.n))/N)*
        ((sum(first(unique(Tabela.n))*(Tabela.Média.-(sum(Tabela.Média)/length(Tabela.n))).^2)/
        (length(Tabela.n)-1)-sum(Tabela.Variância*(first(unique(Tabela.n))-1))/(length(Tabela.n)*
        (first(unique(Tabela.n))-1)))/first(unique(Tabela.n))/length(Tabela.n))+
        ((M-first(unique(Tabela.n)))/M)*(sum(Tabela.Variância*(first(unique(Tabela.n))-1))/(length(Tabela.n)*
        (first(unique(Tabela.n))-1))/
        (length(Tabela.n)*first(unique(Tabela.n))))))),        
        ((sum(Tabela.Média)/length(Tabela.n))+quantile(TDist(length(Tabela.n)-1),1-α/2)*sqrt((((N-length(Tabela.n))/N)*
        ((sum(first(unique(Tabela.n))*(Tabela.Média.-(sum(Tabela.Média)/length(Tabela.n))).^2)/
        (length(Tabela.n)-1)-sum(Tabela.Variância*(first(unique(Tabela.n))-1))/(length(Tabela.n)*
        (first(unique(Tabela.n))-1)))/first(unique(Tabela.n))/length(Tabela.n))+
        ((M-first(unique(Tabela.n)))/M)*(sum(Tabela.Variância*(first(unique(Tabela.n))-1))/(length(Tabela.n)*
        (first(unique(Tabela.n))-1))/(length(Tabela.n)*first(unique(Tabela.n))))))), 
        ((N*M*sum(Tabela.Média)/length(Tabela.n))/Conversor), 
        (((N*M*sum(Tabela.Média)/length(Tabela.n))-(N*M*quantile(TDist(length(Tabela.n)-1),1-α/2)*sqrt((((N-length(Tabela.n))/N)*
        ((sum(first(unique(Tabela.n))*(Tabela.Média.-(sum(Tabela.Média)/length(Tabela.n))).^2)/
        (length(Tabela.n)-1)-sum(Tabela.Variância*(first(unique(Tabela.n))-1))/(length(Tabela.n)*
        (first(unique(Tabela.n))-1)))/first(unique(Tabela.n))/length(Tabela.n))+
        ((M-first(unique(Tabela.n)))/M)*(sum(Tabela.Variância*(first(unique(Tabela.n))-1))/(length(Tabela.n)*
        (first(unique(Tabela.n))-1))/(length(Tabela.n)*first(unique(Tabela.n))))))))/Conversor),
        (((N*M*sum(Tabela.Média)/length(Tabela.n))+(N*M*quantile(TDist(length(Tabela.n)-1),1-α/2)*sqrt((((N-length(Tabela.n))/N)*
        ((sum(first(unique(Tabela.n))*(Tabela.Média.-(sum(Tabela.Média)/length(Tabela.n))).^2)/
        (length(Tabela.n)-1)-sum(Tabela.Variância*(first(unique(Tabela.n))-1))/(length(Tabela.n)*
        (first(unique(Tabela.n))-1)))/first(unique(Tabela.n))/length(Tabela.n))+
        ((M-first(unique(Tabela.n)))/M)*(sum(Tabela.Variância*(first(unique(Tabela.n))-1))/(length(Tabela.n)*
        (first(unique(Tabela.n))-1))/
        (length(Tabela.n)*first(unique(Tabela.n))))))))/Conversor), Area, 
        (quantile(TDist(length(Tabela.n)-1),1-α/2)*sqrt((((N-length(Tabela.n))/N)*((sum(first(unique(Tabela.n))*
        (Tabela.Média.-(sum(Tabela.Média)/length(Tabela.n))).^2)/
        (length(Tabela.n)-1)-sum(Tabela.Variância*(first(unique(Tabela.n))-1))/(length(Tabela.n)*
        (first(unique(Tabela.n))-1)))/first(unique(Tabela.n))/length(Tabela.n))+
        ((M-first(unique(Tabela.n)))/M)*(sum(Tabela.Variância*(first(unique(Tabela.n))-1))/(length(Tabela.n)*
        (first(unique(Tabela.n))-1))/(length(Tabela.n)*first(unique(Tabela.n))))))/
        (sum(Tabela.Média)/length(Tabela.n)))*100, quantile(TDist(length(Tabela.n)-1),1-α/2)*sqrt((((N-length(Tabela.n))/N)*((sum(first(unique(Tabela.n))*
        (Tabela.Média.-(sum(Tabela.Média)/length(Tabela.n))).^2)/
        (length(Tabela.n)-1)-sum(Tabela.Variância*(first(unique(Tabela.n))-1))/
        (length(Tabela.n)*(first(unique(Tabela.n))-1)))/first(unique(Tabela.n))/length(Tabela.n))+
        ((M-first(unique(Tabela.n)))/M)*(sum(Tabela.Variância*(first(unique(Tabela.n))-1))/(length(Tabela.n)*
        (first(unique(Tabela.n))-1))/(length(Tabela.n)*first(unique(Tabela.n)))))), sqrt((((N-length(Tabela.n))/N)*((sum(first(unique(Tabela.n))*(Tabela.Média.-(sum(Tabela.Média)/length(Tabela.n))).^2)/
        (length(Tabela.n)-1)-sum(Tabela.Variância*(first(unique(Tabela.n))-1))/(length(Tabela.n)*
        (first(unique(Tabela.n))-1)))/first(unique(Tabela.n))/length(Tabela.n))+
        ((M-first(unique(Tabela.n)))/M)*(sum(Tabela.Variância*(first(unique(Tabela.n))-1))/(length(Tabela.n)*
        (first(unique(Tabela.n))-1))/(length(Tabela.n)*first(unique(Tabela.n)))))), 
        sum(Tabela.Variância*(first(unique(Tabela.n))-1))/(length(Tabela.n)*(first(unique(Tabela.n))-1))/Conversor, 
        (sum(first(unique(Tabela.n))*(Tabela.Média.-(sum(Tabela.Média)/length(Tabela.n))).^2)/
        (length(Tabela.n)-1)-sum(Tabela.Variância*(first(unique(Tabela.n))-1))/(length(Tabela.n)*
        (first(unique(Tabela.n))-1)))/first(unique(Tabela.n))/Conversor, (((sum(first(unique(Tabela.n))*(Tabela.Média.-(sum(Tabela.Média)/length(Tabela.n))).^2)/
        (length(Tabela.n)-1)-sum(Tabela.Variância*(first(unique(Tabela.n))-1))/(length(Tabela.n)*
        (first(unique(Tabela.n))-1)))/first(unique(Tabela.n))/Conversor)+
        (sum(Tabela.Variância*(first(unique(Tabela.n))-1))/(length(Tabela.n)*(first(unique(Tabela.n))-1))/Conversor)), 
        (((N-length(Tabela.n))/N)*((sum(first(unique(Tabela.n))*(Tabela.Média.-(sum(Tabela.Média)/length(Tabela.n))).^2)/
        (length(Tabela.n)-1)-sum(Tabela.Variância*(first(unique(Tabela.n))-1))/(length(Tabela.n)*
        (first(unique(Tabela.n))-1)))/first(unique(Tabela.n))/length(Tabela.n))+((M-first(unique(Tabela.n)))/M)*
        (sum(Tabela.Variância*(first(unique(Tabela.n))-1))/(length(Tabela.n)*(first(unique(Tabela.n))-1))/
        (length(Tabela.n)*first(unique(Tabela.n)))))/Conversor, (0.1*(sum(Tabela.Média)/length(Tabela.n))), 1-(length(Tabela.n)/N), Tamanho_da_amostra,  
        M, N, length(Tabela.n), first(unique(Tabela.n)), α]), digits = 2)) #Tabela de resultados  

        return Resultados
    end

    function ESTRAT(Dados, Area, AreaParc, EAR, α, estratos, subestrato, selected_variable, unit) #Determinar função
        
        Area = Float64(Meta.parse(Area))
        AreaParc = Float64(Meta.parse(AreaParc))
        α = Float64(Meta.parse(α))
        EAR = Float64(Meta.parse(EAR))

        nh = combine(groupby(Dados, Symbol(estratos)), nrow => :Subestrato_Count)[!,2]

        nh1 = nh[1]
        nh2 = nh[2]
        nh3 = nh[3]

        N = (Area*10000)/AreaParc

        Conversor = 10000/AreaParc
        
        Volume = (Conversor.*Dados[!, selected_variable])
        Estrato = Dados[!,estratos]
        Unidade = Dados[!,subestrato]


        Conjunto_de_dados = DataFrame(Estrato = Estrato, Unidade = Unidade, Volume = Volume)

        Informações_do_inventário = DataFrame(Variáveis=["Área da população (ha)", 
        "Número total potencial de Unidades da população", "Nível de significância (α)", "Número de Unidades amostradas no estrato I", 
        "Número de Unidades amostradas no estrato II", "Número de Unidades amostradas no estrato III", 
        "Número de estratos", "Número de Unidades totais", "Número potencial de Unidades do estrato I", 
        "Número potencial de Unidades do estrato II", "Número potencial de Unidades do estrato III"], 
        Valores=[Area, N, α, nh1, nh2, nh3, length(unique(Estrato)), length(Unidade), 
        (round((Area/(length(Unidade)))*nh1)*10), (round((Area/(length(Unidade)))*nh2)*10), 
        (round((Area/(length(Unidade)))*nh3)*10)]) #Tabela de resultados

    
        #Tabela com estatítica descritiva pro estrato
        Tabela= combine(groupby(Conjunto_de_dados, :Estrato)) do df
            (Unidade=length(unique(df.Unidade)), Total= sum(df.Volume), Média= mean(df.Volume), Variância= var(df.Volume), 
            Erro_padrão= sqrt(var(df.Volume)))
        end

        Anova_da_estratificação = DataFrame(Fontes_de_variação=["Entre estratos", "Dentro dos estratos", "Total"], 
        gl=[length(unique(Estrato))-1 , length(Unidade)-length(unique(Estrato)), length(Unidade)-1], 
        SQ=[sum(Tabela.Unidade.*(Tabela.Média.-mean(Volume)).^2), 
        sum((Volume.-mean(Volume)).^2)-sum(Tabela.Unidade.*(Tabela.Média.-mean(Volume)).^2), 
        sum((Volume.-mean(Volume)).^2)], 
        QM=[(sum(Tabela.Unidade.*(Tabela.Média.-mean(Volume)).^2))/(length(unique(Estrato))-1), 
        (sum((Volume.-mean(Volume)).^2)-sum(Tabela.Unidade.*(Tabela.Média.-mean(Volume)).^2))/
        (length(Unidade)-length(unique(Estrato))), (sum((Volume.-mean(Volume)).^2))/
        (length(Unidade)-1)], F=[((sum(Tabela.Unidade.*(Tabela.Média.-mean(Volume)).^2))/
        (length(unique(Estrato))-1))/
        ((sum((Volume.-mean(Volume)).^2)-sum(Tabela.Unidade.*(Tabela.Média.-mean(Volume)).^2))/
        (length(Unidade)-length(unique(Estrato)))), "", ""])

        t= quantile(TDist(length(Unidade)-1),(1-(α/2))) #Valor de t 

        if (1-(length(Unidade)/N)) ≥ 0.98 #f maior ou igual a 0,98 população infinita
            População = "A população avaliada é considerada infinita"   
        elseif (1-(length(Unidade)/N)) < 0.98 #f menor que 0,98 população finita
            População = "A população avaliada é considerada finita"    
        end

        Tamanho_da_amostra = if (1-(length(Unidade)/N)) ≥ 0.98
            #População infinita. O tamanho da amostra é calculado pela seguinte equação:
            Infinita=(((t)^2)*sum(((round((Area/(length(Unidade)))*nh1)*10)/N; 
            (round((Area/(length(Unidade)))*nh2)*10)/N; 
            (round((Area/(length(Unidade)))*nh3)*10)/N*Tabela.Variância)))/(((0.1*mean(Volume)))^2)
            round(Infinita)
            Finita = Infinita
        elseif (1-(length(Unidade)/N)) < 0.98
            #População finita. O tamanho da amostra é calculado pela seguinte equação:
            Finita=(((t)^2)*sum(((round((Area/(length(Unidade)))*nh1)*10)/N; 
            (round((Area/(length(Unidade)))*nh2)*10)/N; 
            (round((Area/(length(Unidade)))*nh3)*10)/N*Tabela.Variância)))/
            (((0.1*mean(Volume)))^2)+((t)^2)*(sum(((round((Area/(length(Unidade)))*nh1)*10)/N; 
            (round((Area/(length(Unidade)))*nh2)*10)/N; 
            (round((Area/(length(Unidade)))*nh3)*10)/N*Tabela.Variância))/N)
            round(Finita)
        end

        if (((t*sqrt(((((((round((Area/(length(Unidade)))*nh1)*10)/N)^2)*sum(Tabela.Variância/Tabela.Unidade))+
            ((((round((Area/(length(Unidade)))*nh2)*10)/N)^2)*sum(Tabela.Variância/Tabela.Unidade))+
            ((((round((Area/(length(Unidade)))*nh3)*10)/N)^2)*sum(Tabela.Variância/Tabela.Unidade)))/
            length(unique(Estrato)))-(sum(((round((Area/(length(Unidade)))*nh1)*10)/N; (round((Area/(length(Unidade)))*nh2)*10)/N; 
            (round((Area/(length(Unidade)))*nh3)*10)/N*Tabela.Variância))/N)))/mean(Volume))*100) > EAR
            Observação = "Diante do exposto, conclui-se que os resultados obtidos na amostragem não satisfazem as exigências deprecisão estabelecidas para o inventário, ou seja, um erro de amostragem máximo de ±$(EAR)% da média  para confiabilidade designada. \n\nO erro estimado foi maior que o limite fixado, sendo recomendado incluir mais unidades amostrais no inventário."
        elseif (((t*sqrt(((((((round((Area/(length(Unidade)))*nh1)*10)/N)^2)*sum(Tabela.Variância/Tabela.Unidade))+
            ((((round((Area/(length(Unidade)))*nh2)*10)/N)^2)*sum(Tabela.Variância/Tabela.Unidade))+
            ((((round((Area/(length(Unidade)))*nh3)*10)/N)^2)*sum(Tabela.Variância/Tabela.Unidade)))/
            length(unique(Estrato)))-(sum(((round((Area/(length(Unidade)))*nh1)*10)/N; (round((Area/(length(Unidade)))*nh2)*10)/N; 
            (round((Area/(length(Unidade)))*nh3)*10)/N*Tabela.Variância))/N)))/mean(Volume))*100) ≤ EAR
            Observação  = "Diante do exposto, conclui-se que os resultados obtidos na amostragem satisfazem as exigências de precisão estabelecidas para o inventário, ou seja, um erro de amostragem máximo de ±$(EAR)% da média para confiabilidade designada. \n\nO erro estimado foi menor que o limite fixado, assim as unidades amostrais são suficientes para o inventário."
        end  
    

        Resultados = DataFrame(Variáveis=["Média estratificada ($(unit)/ha)", "Limite inferior do intervalo de confiança para Média ($(unit)/ha)", 
        "Limite superior do intervalo de confiança para Média ($(unit)/ha)", "Total da população ($(unit))", 
        "Limite inferior do intervalo de confiança para o total ($(unit))", "Limite superior do intervalo de confiança para o total ($(unit))", "Área da população (ha)",
        "Erro da amostragem relativo (%)", "Erro da amostragem absoluto ($(unit)/ha)", "Erro padrão ($(unit)/ha)", "Desvio padrão ($(unit)/ha)", 
        "Variância estrato I ($(unit)/ha)²", "Variância estrato II ($(unit)/ha)²", "Variância estrato III ($(unit)/ha)²", "Variância estratificada ($(unit)/ha)²", 
        "Variância da média relativa (%)", "Fator de correção", "Limite de erro da amostragem requerido", 
        "Tamanho da amostra estrato I", "Tamanho da amostra estrato II", "Tamanho da amostra estrato III", 
        "Tamanho da amostra", "População", "Observação"], Valores=[mean(Volume), (mean(Volume)-(t*sqrt(((((((round((Area/(length(Unidade)))*nh1)*10)/N)^2)*sum(Tabela.Variância/Tabela.Unidade))+
        ((((round((Area/(length(Unidade)))*nh2)*10)/N)^2)*sum(Tabela.Variância/Tabela.Unidade))+
        ((((round((Area/(length(Unidade)))*nh3)*10)/N)^2)*sum(Tabela.Variância/Tabela.Unidade)))/
        length(unique(Estrato)))-(sum(((round((Area/(length(Unidade)))*nh1)*10)/N; 
        (round((Area/(length(Unidade)))*nh2)*10)/N; 
        (round((Area/(length(Unidade)))*nh3)*10)/
        N*Tabela.Variância))/N)))), (mean(Volume)+(t*sqrt(((((((round((Area/(length(Unidade)))*nh1)*10)/N)^2)*sum(Tabela.Variância/Tabela.Unidade))+
        ((((round((Area/(length(Unidade)))*nh2)*10)/N)^2)*sum(Tabela.Variância/Tabela.Unidade))+
        ((((round((Area/(length(Unidade)))*nh3)*10)/N)^2)*sum(Tabela.Variância/Tabela.Unidade)))/
        length(unique(Estrato)))-(sum(((round((Area/(length(Unidade)))*nh1)*10)/N; 
        (round((Area/(length(Unidade)))*nh2)*10)/N; 
        (round((Area/(length(Unidade)))*nh3)*10)/
        N*Tabela.Variância))/N)))), ((N*mean(Volume))/Conversor), (((N*mean(Volume))-N*(t*sqrt(((((((round((Area/(length(Unidade)))*nh1)*10)/N)^2)*
        sum(Tabela.Variância/Tabela.Unidade))+
        ((((round((Area/(length(Unidade)))*nh2)*10)/N)^2)*sum(Tabela.Variância/Tabela.Unidade))+
        ((((round((Area/(length(Unidade)))*nh3)*10)/N)^2)*sum(Tabela.Variância/Tabela.Unidade)))/
        length(unique(Estrato)))-(sum(((round((Area/(length(Unidade)))*nh1)*10)/N; 
        (round((Area/(length(Unidade)))*nh2)*10)/N; 
        (round((Area/(length(Unidade)))*nh3)*10)/
        N*Tabela.Variância))/N))))/Conversor), (((N*mean(Volume))+N*(t*sqrt(((((((round((Area/(length(Unidade)))*nh1)*10)/N)^2)*
        sum(Tabela.Variância/Tabela.Unidade))+
        ((((round((Area/(length(Unidade)))*nh2)*10)/N)^2)*sum(Tabela.Variância/Tabela.Unidade))+
        ((((round((Area/(length(Unidade)))*nh3)*10)/N)^2)*sum(Tabela.Variância/Tabela.Unidade)))/
        length(unique(Estrato)))-(sum(((round((Area/(length(Unidade)))*nh1)*10)/N; 
        (round((Area/(length(Unidade)))*nh2)*10)/N; 
        (round((Area/(length(Unidade)))*nh3)*10)/
        N*Tabela.Variância))/N))))/Conversor), Area, (((t*sqrt(((((((round((Area/(length(Unidade)))*nh1)*10)/N)^2)*sum(Tabela.Variância/Tabela.Unidade))+
        ((((round((Area/(length(Unidade)))*nh2)*10)/N)^2)*sum(Tabela.Variância/Tabela.Unidade))+
        ((((round((Area/(length(Unidade)))*nh3)*10)/N)^2)*sum(Tabela.Variância/Tabela.Unidade)))/
        length(unique(Estrato)))-(sum(((round((Area/(length(Unidade)))*nh1)*10)/N; 
        (round((Area/(length(Unidade)))*nh2)*10)/N; 
        (round((Area/(length(Unidade)))*nh3)*10)/N*Tabela.Variância))/N)))/
        mean(Volume))*100), (t*sqrt(((((((round((Area/(length(Unidade)))*nh1)*10)/N)^2)*sum(Tabela.Variância/Tabela.Unidade))+
        ((((round((Area/(length(Unidade)))*nh2)*10)/N)^2)*sum(Tabela.Variância/Tabela.Unidade))+
        ((((round((Area/(length(Unidade)))*nh3)*10)/N)^2)*sum(Tabela.Variância/Tabela.Unidade)))/
        length(unique(Estrato)))-(sum(((round((Area/(length(Unidade)))*nh1)*10)/N; 
        (round((Area/(length(Unidade)))*nh2)*10)/N; 
        (round((Area/(length(Unidade)))*nh3)*10)/N*Tabela.Variância))/N))), sqrt(((((((round((Area/(length(Unidade)))*nh1)*10)/N)^2)*sum(Tabela.Variância/Tabela.Unidade))+
        ((((round((Area/(length(Unidade)))*nh2)*10)/N)^2)*sum(Tabela.Variância/Tabela.Unidade))+
        ((((round((Area/(length(Unidade)))*nh3)*10)/N)^2)*sum(Tabela.Variância/Tabela.Unidade)))/
        length(unique(Estrato)))-(sum(((round((Area/(length(Unidade)))*nh1)*10)/N; 
        (round((Area/(length(Unidade)))*nh2)*10)/N; 
        (round((Area/(length(Unidade)))*nh3)*10)/N*Tabela.Variância))/N)), sqrt(mean(Tabela.Variância)), 
        ((((round((Area/(length(Unidade)))*nh1)*10)/N)^2)*sum(Tabela.Variância/Tabela.Unidade)), 
        ((((round((Area/(length(Unidade)))*nh2)*10)/N)^2)*sum(Tabela.Variância/Tabela.Unidade)), 
        ((((round((Area/(length(Unidade)))*nh3)*10)/N)^2)*sum(Tabela.Variância/Tabela.Unidade)), mean(Tabela.Variância), 
        ((((((round((Area/(length(Unidade)))*nh1)*10)/N)^2)*sum(Tabela.Variância/Tabela.Unidade))+
        ((((round((Area/(length(Unidade)))*nh2)*10)/N)^2)*sum(Tabela.Variância/Tabela.Unidade))+
        ((((round((Area/(length(Unidade)))*nh3)*10)/N)^2)*sum(Tabela.Variância/Tabela.Unidade)))/
        length(unique(Estrato)))-(sum(((round((Area/(length(Unidade)))*nh1)*10)/N; 
        (round((Area/(length(Unidade)))*nh2)*10)/N; 
        (round((Area/(length(Unidade)))*nh3)*10)/N*Tabela.Variância))/N), (1-(length(Unidade)/N)), (0.1*mean(Volume)), 
        (round((Area/(length(Unidade)))*nh1)*10)/N*(round(Finita)), 
        (round((Area/(length(Unidade)))*nh2)*10)/N*(round(Finita)), 
        (round((Area/(length(Unidade)))*nh3)*10)/N*(round(Finita)), Tamanho_da_amostra, População, Observação]) #Tabela de resultados    

        #Resultados = [Dados, Informações_do_inventário, Tabela, Anova_da_estratificação, Resultados]"""

        return Resultados

    end#Tabela de resultados   

    function IND(Dados, AreaParc, N1, N2, n_unidade, Ocasiao_1, Ocasiao_2, α, unit)

        AreaParc = Float64(Meta.parse(AreaParc))
        α = Float64(Meta.parse(α))
        N1 = Int64(Meta.parse(N1))
        N2 = Int64(Meta.parse(N2))

        Conversor=1/AreaParc

        ###Primeira ocasião####
        Unidades = Dados[!, n_unidade]
        Ocasião_1 = (Conversor.*Dados[!, Ocasiao_1])
        Ocasião_2 = (Conversor.*Dados[!, Ocasiao_2])
        Independente = DataFrame(Unidades = Dados.n, Ocasião_1 = Ocasião_1, Ocasião_2 = Ocasião_2)
        
        Primeira_ocasião = DataFrame(Variáveis=["Média ($(unit)/ha)", "Limite inferior do intervalo de confiança para média ($(unit)/ha)", 
        "Limite superior do intervalo de confiança para média ($(unit)/ha)", "Total da população ($(unit))", "Limite inferior do intervalo de confiança para o total ($(unit))", 
        "Limite superior do intervalo de confiança para o total ($(unit))", "Área da população (ha)", "Erro da amostragem relativo (%)", 
        "Erro da amostragem absoluto ($(unit)/ha)", "Erro padrão ($(unit)/ha)", "Desvio padrão ($(unit)/ha)", "Variância ($(unit)/ha)²", "Variância da média ($(unit)/ha)²", 
        "Número total de unidades", "Nível de de significância (α)"], 
        Valores=[ mean(Ocasião_1),  (mean(Ocasião_1)-quantile(TDist(length(Ocasião_1)-1),1-α/2)*sqrt.(((sum((Ocasião_1.-mean(Ocasião_1)).^2)/
        length(Ocasião_1)-1)/length(Ocasião_1))*(1-(length(Ocasião_1)/N1)))), (mean(Ocasião_1)+quantile(TDist(length(Ocasião_1)-1),1-α/2)*
        sqrt.(((sum((Ocasião_1.-mean(Ocasião_1)).^2)/length(Ocasião_1)-1)/length(Ocasião_1))*(1-(length(Ocasião_1)/N1)))), ((N1*mean(Ocasião_1))/Conversor),
        (((N1*mean(Ocasião_1))-N1*quantile(TDist(length(Ocasião_1)-1),1-α/2)*sqrt.(((sum((Ocasião_1.-mean(Ocasião_1)).^2)/
        length(Ocasião_1)-1)/length(Ocasião_1))* (1-(length(Ocasião_1)/N1))))/Conversor), (((N1*mean(Ocasião_1))+N1*quantile(TDist(length(Ocasião_1)-1),1-α/2)*
        sqrt.(((sum((Ocasião_1.-mean(Ocasião_1)).^2)/length(Ocasião_1)-1)/length(Ocasião_1))*(1-(length(Ocasião_1)/N1))))/Conversor), N1, 
        (quantile(TDist(length(Ocasião_1)-1),1-α/2)*sqrt.(((sum((Ocasião_1.-mean(Ocasião_1)).^2)/length(Ocasião_1)-1)/
        length(Ocasião_1))*(1-(length(Ocasião_1)/N1)))/mean(Ocasião_1))*100, quantile(TDist(length(Ocasião_1)-1),1-α/2)*sqrt.(((sum((Ocasião_1.-mean(Ocasião_1)).^2)/length(Ocasião_1)-1)/
        length(Ocasião_1))*(1-(length(Ocasião_1)/N1))), sqrt.(((sum((Ocasião_1.-mean(Ocasião_1)).^2)/length(Ocasião_1)-1)/length(Ocasião_1))*
        (1-(length(Ocasião_1)/N1))), sqrt(sum((Ocasião_1.-mean(Ocasião_1)).^2)/length(Ocasião_1)-1), sum((Ocasião_1.-mean(Ocasião_1)).^2)/length(Ocasião_1)-1, 
        ((sum((Ocasião_1.-mean(Ocasião_1)).^2)/length(Ocasião_1)-1)/length(Ocasião_1))*
        (1-(length(Ocasião_1)/N1)), length(Ocasião_1), α]) #Tabela de resultados 
        
        Segunda_ocasião = DataFrame(Variáveis=["Média ($(unit)/ha)", "Limite inferior do intervalo de confiança para média ($(unit)/ha)", 
        "Limite superior do intervalo de confiança para média ($(unit)/ha)", "Total da população ($(unit))", "Limite inferior do intervalo de confiança para o total ($(unit))", 
        "Limite superior do intervalo de confiança para o total ($(unit))", "Área da população (ha)", "Erro da amostragem relativo (%)", 
        "Erro da amostragem absoluto ($(unit)/ha)", "Erro padrão ($(unit)/ha)", "Desvio padrão ($(unit)/ha)", "Variância ($(unit)/ha)²", "Variância da média ($(unit)/ha)²", 
        "Número total de unidades", "Nível de de significância (α)"], Valores=[mean(Ocasião_2), ((mean(Ocasião_2))-(quantile(TDist((length(Ocasião_2))-1),1-α/2))*
        sqrt.(((sum((Ocasião_2.-(mean(Ocasião_2))).^2)/length(Ocasião_2)-1)/(length(Ocasião_2)))*(1-((length(Ocasião_2))/N2)))), 
        ((mean(Ocasião_2))+(quantile(TDist((length(Ocasião_2))-1),1-α/2))*sqrt.(((sum((Ocasião_2.-(mean(Ocasião_2))).^2)/
        length(Ocasião_2)-1)/(length(Ocasião_2)))*(1-((length(Ocasião_2))/N2)))), ((N2*mean(Ocasião_2))/Conversor), 
        ((N2*(mean(Ocasião_2)))-N2*((quantile(TDist((length(Ocasião_2))-1),1-α/2))*sqrt.(((sum((Ocasião_2.-(mean(Ocasião_2))).^2)/
        length(Ocasião_2)-1)/(length(Ocasião_2)))*(1-((length(Ocasião_2))/N2))))/Conversor), (((N2*(mean(Ocasião_2)))+
        N2*(quantile(TDist((length(Ocasião_2))-1),1-α/2))*sqrt.(((sum((Ocasião_2.-(mean(Ocasião_2))).^2)/
        length(Ocasião_2)-1)/(length(Ocasião_2)))*(1-((length(Ocasião_2))/N2))))/Conversor), N2, 
        ((quantile(TDist((length(Ocasião_2))-1),1-α/2))*sqrt.(((sum((Ocasião_2.-(mean(Ocasião_2))).^2)/
        length(Ocasião_2)-1)/(length(Ocasião_2)))*(1-((length(Ocasião_2))/N2)))/(mean(Ocasião_2)))*100, 
        (quantile(TDist((length(Ocasião_2))-1),1-α/2))*sqrt.(((sum((Ocasião_2.-(mean(Ocasião_2))).^2)/
        length(Ocasião_2)-1)/(length(Ocasião_2)))*(1-((length(Ocasião_2))/N2))), sqrt.(((sum((Ocasião_2.-(mean(Ocasião_2))).^2)/
        length(Ocasião_2)-1)/(length(Ocasião_2)))*(1-((length(Ocasião_2))/N2))), sqrt(sum((Ocasião_2.-(mean(Ocasião_2))).^2)/length(Ocasião_2)-1), 
        (sum((Ocasião_2.-(mean(Ocasião_2))).^2)/length(Ocasião_2)-1), ((sum((Ocasião_2.-(mean(Ocasião_2))).^2)/length(Ocasião_2)-1)/
        (length(Ocasião_2)))*(1-((length(Ocasião_2))/N2)), length(Ocasião_2), α]) #Tabela de resultados
    
        
        Mudança_crescimento = DataFrame(Variáveis=["Crescimento médio ($(unit)/ha)", "Limite inferior do intervalo de confiança para média ($(unit)/ha)", 
        "Limite superior do intervalo de confiança para média ($(unit)/ha)", "Crescimento total estimado ($(unit))", "Limite inferior do intervalo de confiança para o total ($(unit))", 
        "Limite superior do intervalo de confiança para o total ($(unit))", "Área da população (ha)", "Erro da amostragem relativo (%)", 
        "Erro da amostragem absoluto ($(unit)/ha)", "Erro padrão ($(unit)/ha)", "Variância da média ($(unit)/ha)²"], Valores=[mean(Ocasião_2)-mean(Ocasião_1), 
        (mean(Ocasião_2)-mean(Ocasião_1)-((quantile(TDist((length(Ocasião_1)-1)+((length(Ocasião_2))-1)-1),1-α/2))*
        (sqrt((((sum((Ocasião_1.-mean(Ocasião_1)).^2)/length(Ocasião_1)-1)/length(Ocasião_1))*
        (1-(length(Ocasião_1)/N1)))+(((sum((Ocasião_2.-(mean(Ocasião_2))).^2)/length(Ocasião_2)-1)/
        (length(Ocasião_2)))*(1-((length(Ocasião_2))/N2))))))), (mean(Ocasião_2)-mean(Ocasião_1)+((quantile(TDist((length(Ocasião_1)-1)+
        ((length(Ocasião_2))-1)-1),1-α/2))*(sqrt((((sum((Ocasião_1.-mean(Ocasião_1)).^2)/length(Ocasião_1)-1)/length(Ocasião_1))*
        (1-(length(Ocasião_1)/N1)))+(((sum((Ocasião_2.-(mean(Ocasião_2))).^2)/length(Ocasião_2)-1)/
        (length(Ocasião_2)))*(1-((length(Ocasião_2))/N2))))))), N2*(mean(Ocasião_2)-mean(Ocasião_1)), 
        ((N2*(mean(Ocasião_2)-mean(Ocasião_1)))-N2*((quantile(TDist((length(Ocasião_1)-1)+
        ((length(Ocasião_2))-1)-1),1-α/2))*(sqrt((((sum((Ocasião_1.-mean(Ocasião_1)).^2)/
        length(Ocasião_1)-1)/length(Ocasião_1))*(1-(length(Ocasião_1)/N1)))+
        (((sum((Ocasião_2.-(mean(Ocasião_2))).^2)/length(Ocasião_2)-1)/(length(Ocasião_2)))*
        (1-((length(Ocasião_2))/N2))))))), ((N2*(mean(Ocasião_2)-mean(Ocasião_1)))+N2*((quantile(TDist((length(Ocasião_1)-1)+
        ((length(Ocasião_2))-1)-1),1-α/2))*(sqrt((((sum((Ocasião_1.-mean(Ocasião_1)).^2)/
        length(Ocasião_1)-1)/length(Ocasião_1))*(1-(length(Ocasião_1)/N1)))+(((sum((Ocasião_2.-(mean(Ocasião_2))).^2)/length(Ocasião_2)-1)/
        (length(Ocasião_2)))*(1-((length(Ocasião_2))/N2))))))), N2, (((quantile(TDist((length(Ocasião_1)-1)+((length(Ocasião_2))-1)-1),1-α/2))*
        (sqrt((((sum((Ocasião_1.-mean(Ocasião_1)).^2)/length(Ocasião_1)-1)/length(Ocasião_1))*
        (1-(length(Ocasião_1)/N1)))+(((sum((Ocasião_2.-(mean(Ocasião_2))).^2)/length(Ocasião_2)-1)/
        (length(Ocasião_2)))*(1-((length(Ocasião_2))/N2))))))/(mean(Ocasião_2)-mean(Ocasião_1)))*100, 
        ((quantile(TDist((length(Ocasião_1)-1)+((length(Ocasião_2))-1)-1),1-α/2))*
        (sqrt((((sum((Ocasião_1.-mean(Ocasião_1)).^2)/length(Ocasião_1)-1)/length(Ocasião_1))*
        (1-(length(Ocasião_1)/N1)))+(((sum((Ocasião_2.-(mean(Ocasião_2))).^2)/length(Ocasião_2)-1)/
        (length(Ocasião_2)))*(1-((length(Ocasião_2))/N2)))))), sqrt((((sum((Ocasião_1.-mean(Ocasião_1)).^2)/length(Ocasião_1)-1)/length(Ocasião_1))*
        (1-(length(Ocasião_1)/N1)))+(((sum((Ocasião_2.-(mean(Ocasião_2))).^2)/length(Ocasião_2)-1)/(length(Ocasião_2)))*
        (1-((length(Ocasião_2))/N2)))), (((sum((Ocasião_1.-mean(Ocasião_1)).^2)/length(Ocasião_1)-1)/length(Ocasião_1))*(1-(length(Ocasião_1)/N1)))+
        (((sum((Ocasião_2.-(mean(Ocasião_2))).^2)/length(Ocasião_2)-1)/(length(Ocasião_2)))*(1-((length(Ocasião_2))/N2)))]) #Tabela de resultados  

        Resultados = [Independente, Primeira_ocasião, Segunda_ocasião, Mudança_crescimento]
        
        return Mudança_crescimento

    end

    function MULTI(Dados, AreaParc, Area, EAR, α, unit) #Determinar função

        Area = Float64(Meta.parse(Area))
        AreaParc = Float64(Meta.parse(AreaParc))
        α = Float64(Meta.parse(α))
        EAR = Float64(Meta.parse(EAR))
        
        N = (Area)/2

        Conversor = 10000/AreaParc

        Conjunto_de_dados = (Conversor.*Dados)
        
        #Tabela com estatísticas descritivas por conglomerados
        Tabela=transform(Conjunto_de_dados, AsTable(:) .=> ByRow.([I -> count(!ismissing, I), sum, mean, var]).=>[:n, :Soma, :Média, :Variância])
        
        if  (((quantile(TDist((length(Tabela.n))-1),1-α/2))*(sqrt.((((sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/((length(Tabela.n))*
            (((first(unique((Tabela.n))))).-1)))+(((sum(((first(unique((Tabela.n))))).*(Tabela.Média.-(sum(Tabela.Média)/
            (length(Tabela.n)))).^2)/((length(Tabela.n))-1)).-
            (sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/((length(Tabela.n))*
            (((first(unique((Tabela.n))))).-1))))./((first(unique((Tabela.n)))))))./((length(Tabela.n))*
            ((first(unique((Tabela.n)))))))*(1 .+((((sum(((first(unique((Tabela.n))))).*(Tabela.Média.-(sum(Tabela.Média)/
            (length(Tabela.n)))).^2)/((length(Tabela.n))-1)).-
            (sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/((length(Tabela.n))*
            (((first(unique((Tabela.n))))).-1))))./((first(unique((Tabela.n))))))/
            ((((sum(((first(unique((Tabela.n))))).*(Tabela.Média.-(sum(Tabela.Média)/(length(Tabela.n)))).^2)/
            ((length(Tabela.n))-1)).-(sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/
            ((length(Tabela.n))*(((first(unique((Tabela.n))))).-1))))./((first(unique((Tabela.n))))))+
            (sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/((length(Tabela.n))*(((first(unique((Tabela.n))))).-1)))))*
            (((first(unique((Tabela.n))))).-1)))))/(sum(Tabela.Média)/(length(Tabela.n))))*100 > EAR
            Observação = "Diante do exposto, conclui-se que os resultados obtidos na amostragem não satisfazem as exigências de precisão estabelecidas para o inventário, ou seja, um erro de amostragem máximo de $(Int(EAR))% da média  para confiabilidade designada. \n\nO erro estimado foi maior que o limite fixado, sendo recomendado incluir mais unidades amostrais no inventário."
        elseif  (((quantile(TDist((length(Tabela.n))-1),1-α/2))*(sqrt.((((sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/((length(Tabela.n))*
            (((first(unique((Tabela.n))))).-1)))+(((sum(((first(unique((Tabela.n))))).*(Tabela.Média.-(sum(Tabela.Média)/
            (length(Tabela.n)))).^2)/((length(Tabela.n))-1)).-
            (sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/((length(Tabela.n))*
            (((first(unique((Tabela.n))))).-1))))./((first(unique((Tabela.n)))))))./((length(Tabela.n))*
            ((first(unique((Tabela.n)))))))*(1 .+((((sum(((first(unique((Tabela.n))))).*(Tabela.Média.-(sum(Tabela.Média)/
            (length(Tabela.n)))).^2)/((length(Tabela.n))-1)).-
            (sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/((length(Tabela.n))*
            (((first(unique((Tabela.n))))).-1))))./((first(unique((Tabela.n))))))/
            ((((sum(((first(unique((Tabela.n))))).*(Tabela.Média.-(sum(Tabela.Média)/(length(Tabela.n)))).^2)/
            ((length(Tabela.n))-1)).-(sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/
            ((length(Tabela.n))*(((first(unique((Tabela.n))))).-1))))./((first(unique((Tabela.n))))))+
            (sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/((length(Tabela.n))*(((first(unique((Tabela.n))))).-1)))))*
            (((first(unique((Tabela.n))))).-1)))))/(sum(Tabela.Média)/(length(Tabela.n))))*100 ≤ EAR
            Observação  = "Diante do exposto, conclui-se que os resultados obtidos na amostragem satisfazem as exigências de precisão estabelecidas para o inventário, ou seja, um erro de amostragem máximo de ±$(Int(EAR))% da média para confiabilidade designada. \n\nO erro estimado foi menor que o limite fixado, assim as unidades amostrais são suficientes para o inventário."
        end

        Resultados = DataFrame(Variáveis=["Média ($(unit)/ha)", "Limite inferior do intervalo de confiança para média ($(unit)/ha)", 
        "Limite superior do intervalo de confiança para média ($(unit)/ha)", "Total da população ($(unit))", "Limite inferior do intervalo de confiança para o total ($(unit))", 
        "Limite superior do intervalo de confiança para o total ($(unit))", "Área da população (ha)", "Erro da amostragem relativo (%)", "Erro da amostragem absoluto ($(unit)/ha)", 
        "Erro padrão ($(unit)/ha)", "Variância dentro dos conglomerados ($(unit)/ha)²", "Variância entre conglomerados ($(unit)/ha)²", "Variância total ($(unit)/ha)²", "Variância da média ($(unit)/ha)²", 
        "Coeficiente de correlação intraconglomerados", "Tamanho da amostra", "Limite do erro de amostragem requerido", "Número de conglomerados", 
        "Número de subunidades", "Nível de significância (α)", "Observação"], 
        Valores=[sum(Tabela.Média)/(length(Tabela.n)), ((sum(Tabela.Média)/(length(Tabela.n))).-((quantile(TDist((length(Tabela.n))-1),1-α/2))*
        (sqrt.((((sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/((length(Tabela.n))*(((first(unique((Tabela.n))))).-1)))+
        (((sum(((first(unique((Tabela.n))))).*(Tabela.Média.-(sum(Tabela.Média)/(length(Tabela.n)))).^2)/
        ((length(Tabela.n))-1)).-(sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/
        ((length(Tabela.n))*(((first(unique((Tabela.n))))).-1))))./((first(unique((Tabela.n)))))))./
        ((length(Tabela.n))*((first(unique((Tabela.n)))))))*(1 .+((((sum(((first(unique((Tabela.n))))).*
        (Tabela.Média.-(sum(Tabela.Média)/(length(Tabela.n)))).^2)/
        ((length(Tabela.n))-1)).-(sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/
        ((length(Tabela.n))*(((first(unique((Tabela.n))))).-1))))./((first(unique((Tabela.n))))))/
        ((((sum(((first(unique((Tabela.n))))).*(Tabela.Média.-(sum(Tabela.Média)/(length(Tabela.n)))).^2)/
        ((length(Tabela.n))-1)).-(sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/
        ((length(Tabela.n))*(((first(unique((Tabela.n))))).-1))))./((first(unique((Tabela.n))))))+
        (sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/((length(Tabela.n))*(((first(unique((Tabela.n))))).-1)))))*
        (((first(unique((Tabela.n))))).-1)))))), ((sum(Tabela.Média)/(length(Tabela.n))).+((quantile(TDist((length(Tabela.n))-1),1-α/2))*
        (sqrt.((((sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/((length(Tabela.n))*(((first(unique((Tabela.n))))).-1)))+
        (((sum(((first(unique((Tabela.n))))).*(Tabela.Média.-(sum(Tabela.Média)/(length(Tabela.n)))).^2)/
        ((length(Tabela.n))-1)).-(sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/
        ((length(Tabela.n))*(((first(unique((Tabela.n))))).-1))))./((first(unique((Tabela.n)))))))./
        ((length(Tabela.n))*((first(unique((Tabela.n)))))))*(1 .+((((sum(((first(unique((Tabela.n))))).*
        (Tabela.Média.-(sum(Tabela.Média)/(length(Tabela.n)))).^2)/
        ((length(Tabela.n))-1)).-(sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/
        ((length(Tabela.n))*(((first(unique((Tabela.n))))).-1))))./((first(unique((Tabela.n))))))/
        ((((sum(((first(unique((Tabela.n))))).*(Tabela.Média.-(sum(Tabela.Média)/(length(Tabela.n)))).^2)/
        ((length(Tabela.n))-1)).-(sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/
        ((length(Tabela.n))*(((first(unique((Tabela.n))))).-1))))./
        ((first(unique((Tabela.n))))))+(sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/((length(Tabela.n))*
        (((first(unique((Tabela.n))))).-1)))))*
        (((first(unique((Tabela.n))))).-1)))))), ((N*((first(unique((Tabela.n)))))*(sum(Tabela.Média)/(length(Tabela.n))))/Conversor), 
        (((N*((first(unique((Tabela.n)))))*(sum(Tabela.Média)/(length(Tabela.n))))-(N*((first(unique((Tabela.n))))).*
        ((quantile(TDist((length(Tabela.n))-1),1-α/2))*(sqrt.((((sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/
        ((length(Tabela.n))*(((first(unique((Tabela.n))))).-1)))+
        (((sum(((first(unique((Tabela.n))))).*(Tabela.Média.-(sum(Tabela.Média)/(length(Tabela.n)))).^2)/
        ((length(Tabela.n))-1)).-(sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/
        ((length(Tabela.n))*(((first(unique((Tabela.n))))).-1))))./((first(unique((Tabela.n)))))))./
        ((length(Tabela.n))*((first(unique((Tabela.n)))))))*(1 .+((((sum(((first(unique((Tabela.n))))).*
        (Tabela.Média.-(sum(Tabela.Média)/(length(Tabela.n)))).^2)/
        ((length(Tabela.n))-1)).-(sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/
        ((length(Tabela.n))*(((first(unique((Tabela.n))))).-1))))./((first(unique((Tabela.n))))))/
        ((((sum(((first(unique((Tabela.n))))).*(Tabela.Média.-(sum(Tabela.Média)/(length(Tabela.n)))).^2)/
        ((length(Tabela.n))-1)).-(sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/
        ((length(Tabela.n))*(((first(unique((Tabela.n))))).-1))))./((first(unique((Tabela.n))))))+
        (sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/((length(Tabela.n))*
        (((first(unique((Tabela.n))))).-1)))))*
        (((first(unique((Tabela.n))))).-1)))))))/Conversor), (((N*((first(unique((Tabela.n)))))*(sum(Tabela.Média)/(length(Tabela.n))))+(N*((first(unique((Tabela.n))))).*
        ((quantile(TDist((length(Tabela.n))-1),1-α/2))*(sqrt.((((sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/
        ((length(Tabela.n))*(((first(unique((Tabela.n))))).-1)))+
        (((sum(((first(unique((Tabela.n))))).*(Tabela.Média.-(sum(Tabela.Média)/(length(Tabela.n)))).^2)/
        ((length(Tabela.n))-1)).-(sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/
        ((length(Tabela.n))*(((first(unique((Tabela.n))))).-1))))./((first(unique((Tabela.n)))))))./
        ((length(Tabela.n))*((first(unique((Tabela.n)))))))*(1 .+((((sum(((first(unique((Tabela.n))))).*
        (Tabela.Média.-(sum(Tabela.Média)/(length(Tabela.n)))).^2)/
        ((length(Tabela.n))-1)).-(sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/
        ((length(Tabela.n))*(((first(unique((Tabela.n))))).-1))))./((first(unique((Tabela.n))))))/
        ((((sum(((first(unique((Tabela.n))))).*(Tabela.Média.-(sum(Tabela.Média)/(length(Tabela.n)))).^2)/
        ((length(Tabela.n))-1)).-(sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/
        ((length(Tabela.n))*(((first(unique((Tabela.n))))).-1))))./((first(unique((Tabela.n))))))+
        (sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/((length(Tabela.n))*
        (((first(unique((Tabela.n))))).-1)))))*(((first(unique((Tabela.n))))).-1)))))))/Conversor), 
        Area, (((quantile(TDist((length(Tabela.n))-1),1-α/2))*(sqrt.((((sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/
        ((length(Tabela.n))*(((first(unique((Tabela.n))))).-1)))+(((sum(((first(unique((Tabela.n))))).*(Tabela.Média.-(sum(Tabela.Média)/
        (length(Tabela.n)))).^2)/((length(Tabela.n))-1)).-(sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/((length(Tabela.n))*
        (((first(unique((Tabela.n))))).-1))))./((first(unique((Tabela.n)))))))./((length(Tabela.n))*
        ((first(unique((Tabela.n)))))))*(1 .+((((sum(((first(unique((Tabela.n))))).*(Tabela.Média.-(sum(Tabela.Média)/
        (length(Tabela.n)))).^2)/((length(Tabela.n))-1)).-(sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/((length(Tabela.n))*
        (((first(unique((Tabela.n))))).-1))))./((first(unique((Tabela.n))))))/
        ((((sum(((first(unique((Tabela.n))))).*(Tabela.Média.-(sum(Tabela.Média)/(length(Tabela.n)))).^2)/
        ((length(Tabela.n))-1)).-(sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/
        ((length(Tabela.n))*(((first(unique((Tabela.n))))).-1))))./((first(unique((Tabela.n))))))+
        (sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/((length(Tabela.n))*(((first(unique((Tabela.n))))).-1)))))*
        (((first(unique((Tabela.n))))).-1)))))/(sum(Tabela.Média)/(length(Tabela.n))))*100, (quantile(TDist((length(Tabela.n))-1),1-α/2))*(sqrt.((((sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/((length(Tabela.n))*
        (((first(unique((Tabela.n))))).-1)))+(((sum(((first(unique((Tabela.n))))).*(Tabela.Média.-(sum(Tabela.Média)/
        (length(Tabela.n)))).^2)/((length(Tabela.n))-1)).-(sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/((length(Tabela.n))*
        (((first(unique((Tabela.n))))).-1))))./((first(unique((Tabela.n)))))))./((length(Tabela.n))*
        ((first(unique((Tabela.n)))))))*(1 .+((((sum(((first(unique((Tabela.n))))).*
        (Tabela.Média.-(sum(Tabela.Média)/(length(Tabela.n)))).^2)/((length(Tabela.n))-1)).-
        (sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/((length(Tabela.n))*(((first(unique((Tabela.n))))).-1))))./
        ((first(unique((Tabela.n))))))/((((sum(((first(unique((Tabela.n))))).*(Tabela.Média.-(sum(Tabela.Média)/(length(Tabela.n)))).^2)/
        ((length(Tabela.n))-1)).-(sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/((length(Tabela.n))*(((first(unique((Tabela.n))))).-1))))./
        ((first(unique((Tabela.n))))))+(sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/((length(Tabela.n))*
        (((first(unique((Tabela.n))))).-1)))))*(((first(unique((Tabela.n))))).-1)))), sqrt.((((sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/((length(Tabela.n))*(((first(unique((Tabela.n))))).-1)))+
        (((sum(((first(unique((Tabela.n))))).*(Tabela.Média.-(sum(Tabela.Média)/(length(Tabela.n)))).^2)/
        ((length(Tabela.n))-1)).-(sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/((length(Tabela.n))*
        (((first(unique((Tabela.n))))).-1))))./((first(unique((Tabela.n)))))))./((length(Tabela.n))*
        ((first(unique((Tabela.n)))))))*(1 .+((((sum(((first(unique((Tabela.n))))).*(Tabela.Média.-(sum(Tabela.Média)/
        (length(Tabela.n)))).^2)/((length(Tabela.n))-1)).-(sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/((length(Tabela.n))*
        (((first(unique((Tabela.n))))).-1))))./((first(unique((Tabela.n))))))/((((sum(((first(unique((Tabela.n))))).*
        (Tabela.Média.-(sum(Tabela.Média)/(length(Tabela.n)))).^2)/
        ((length(Tabela.n))-1)).-(sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/((length(Tabela.n))*
        (((first(unique((Tabela.n))))).-1))))./((first(unique((Tabela.n))))))+
        (sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/((length(Tabela.n))*(((first(unique((Tabela.n))))).-1)))))*
        (((first(unique((Tabela.n))))).-1))), sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/((length(Tabela.n))*(((first(unique((Tabela.n))))).-1))/Conversor, 
        ((sum(((first(unique((Tabela.n))))).*(Tabela.Média.-(sum(Tabela.Média)/(length(Tabela.n)))).^2)/
        ((length(Tabela.n))-1)).-(sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/((length(Tabela.n))*
        (((first(unique((Tabela.n))))).-1))))./((first(unique((Tabela.n)))))/Conversor, ((sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/((length(Tabela.n))*(((first(unique((Tabela.n))))).-1)))/Conversor)+
        ((((sum(((first(unique((Tabela.n))))).*(Tabela.Média.-(sum(Tabela.Média)/(length(Tabela.n)))).^2)/
        ((length(Tabela.n))-1)).-(sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/((length(Tabela.n))*
        (((first(unique((Tabela.n))))).-1))))./((first(unique((Tabela.n))))))/Conversor), ((((sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/((length(Tabela.n))*(((first(unique((Tabela.n))))).-1)))/Conversor)+
        ((((sum(((first(unique((Tabela.n))))).*(Tabela.Média.-(sum(Tabela.Média)/(length(Tabela.n)))).^2)/
        ((length(Tabela.n))-1)).-(sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/((length(Tabela.n))*
        (((first(unique((Tabela.n))))).-1))))./((first(unique((Tabela.n))))))/Conversor))./((length(Tabela.n))*
        ((first(unique((Tabela.n)))))))*(1 .+((((sum(((first(unique((Tabela.n))))).*(Tabela.Média.-(sum(Tabela.Média)/
        (length(Tabela.n)))).^2)/((length(Tabela.n))-1)).-(sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/((length(Tabela.n))*
        (((first(unique((Tabela.n))))).-1))))./((first(unique((Tabela.n))))))/
        ((((sum(((first(unique((Tabela.n))))).*(Tabela.Média.-(sum(Tabela.Média)/(length(Tabela.n)))).^2)/
        ((length(Tabela.n))-1)).-(sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/((length(Tabela.n))*
        (((first(unique((Tabela.n))))).-1))))./((first(unique((Tabela.n))))))+
        (sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/((length(Tabela.n))*(((first(unique((Tabela.n))))).-1)))))*
        (((first(unique((Tabela.n))))).-1)), (((sum(((first(unique((Tabela.n))))).*(Tabela.Média.-(sum(Tabela.Média)/(length(Tabela.n)))).^2)/
        ((length(Tabela.n))-1)).-(sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/((length(Tabela.n))*
        (((first(unique((Tabela.n))))).-1))))./((first(unique((Tabela.n))))))/
        ((((sum(((first(unique((Tabela.n))))).*(Tabela.Média.-(sum(Tabela.Média)/(length(Tabela.n)))).^2)/
        ((length(Tabela.n))-1)).-(sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/((length(Tabela.n))*
        (((first(unique((Tabela.n))))).-1))))./((first(unique((Tabela.n))))))+(sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/
        ((length(Tabela.n))*(((first(unique((Tabela.n))))).-1)))), (((((quantile(TDist((length(Tabela.n))-1),1-α/2))^2)*((sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/
        ((length(Tabela.n))*(((first(unique((Tabela.n))))).-1)))+
        (((sum(((first(unique((Tabela.n))))).*(Tabela.Média.-(sum(Tabela.Média)/(length(Tabela.n)))).^2)/
        ((length(Tabela.n))-1)).-(sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/((length(Tabela.n))*
        (((first(unique((Tabela.n))))).-1))))./((first(unique((Tabela.n)))))))))./
        ((((0.1*(sum(Tabela.Média)/(length(Tabela.n))))).^2).*((first(unique((Tabela.n))))))).*(1 .+((((sum(((first(unique((Tabela.n))))).*
        (Tabela.Média.-(sum(Tabela.Média)/(length(Tabela.n)))).^2)/
        ((length(Tabela.n))-1)).-(sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/((length(Tabela.n))*
        (((first(unique((Tabela.n))))).-1))))./((first(unique((Tabela.n))))))/
        ((((sum(((first(unique((Tabela.n))))).*(Tabela.Média.-(sum(Tabela.Média)/(length(Tabela.n)))).^2)/
        ((length(Tabela.n))-1)).-(sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/((length(Tabela.n))*
        (((first(unique((Tabela.n))))).-1))))./
        ((first(unique((Tabela.n))))))+(sum(Tabela.Variância.*(((first(unique((Tabela.n))))).-1))/((length(Tabela.n))*
        (((first(unique((Tabela.n))))).-1)))))*(((first(unique((Tabela.n))))).-1)), (0.1*(sum(Tabela.Média)/(length(Tabela.n)))), 
        length(Tabela.n), first((first(unique((Tabela.n))))), α, Observação]) #Tabela de resultados    
        
        #Resultados = [Conjunto_de_dados, Tabela, Resultados]

        return Conjunto_de_dados

    end

    function SIST(Dados, AreaParc, Area, EAR, α, unit) #Determinar função

        Area = Float64(Meta.parse(Area))
        AreaParc = Float64(Meta.parse(AreaParc))
        α = Float64(Meta.parse(α))
        EAR = Float64(Meta.parse(EAR))

        N = (Area*10000)/AreaParc

        Conversor = 10000/AreaParc

        Conjunto_de_dados = (Conversor.*Dados)
        #Tabela com estatítica descritiva por unidade secundária/bloco
        Tabela=transform(Conjunto_de_dados, AsTable(:) .=> ByRow.([I -> count(!ismissing, I), sum, mean, var]).=>[:n, :Soma, :Média, :Variância])
        
        if (1-((length(Tabela.n)*first(unique((Tabela.n))))/N)) ≥ 0.98 #f maior ou igual a 0,98 população infinita
            População = "A população avalaida é considerada infinita"   
        elseif (1-((length(Tabela.n)*first(unique((Tabela.n))))/N)) < 0.98 #f menor que 0,98 população finita
            População = "A população avaliada é considerada finita"    
        end
        
        g=Matrix(Dados)
        matriz=transpose(g)
        
        global a = 0
        global Sx² = 0

        for i in 1:length(matriz)-1
            if i % first(unique((Tabela.n))) == 0
                continue
            end

            global a += (sum(matriz[i]-matriz[i+1])^2)

            Sx² = (a/(2*(length(Tabela.n)*first(unique((Tabela.n))))*((length(Tabela.n)*first(unique((Tabela.n))))-length(Tabela.n))))*
                (1-((length(Tabela.n)*first(unique((Tabela.n))))/N)) #Estimativa aproximada da variância da média
        end

        # Verifica se a população é FINITA ou INFINITA
        if (((quantile(TDist((length(Tabela.n)*first(unique((Tabela.n))))-1), 1-α/2))*sqrt(Sx²))/mean(Tabela.Média))*100 > EAR
            Observação = "Diante do exposto, conclui-se que os resultados obtidos na amostragem não satisfazem as exigências de precisão estabelecidas para o inventário, ou seja, um erro de amostragem máximo de $(Int(EAR))% da média  para confiabilidade designada. \n\nO erro estimado foi maior que o limite fixado, sendo recomendado incluir mais unidades amostrais no inventário."
        elseif (((quantile(TDist((length(Tabela.n)*first(unique((Tabela.n))))-1), 1-α/2))*sqrt(Sx²))/mean(Tabela.Média))*100 ≤ EAR
            Observação  = "Diante do exposto, conclui-se que os resultados obtidos na amostragem satisfazem as exigências de precisão estabelecidas para o inventário, ou seja, um erro de amostragem máximo de ±$(Int(EAR))% da média para confiabilidade designada. \n\nO erro estimado foi menor que o limite fixado, assim as unidades amostrais são suficientes para o inventário."
        end

        Resultados = (
            DataFrame(Variáveis=["Média ($(unit)/ha)", "Limite inferior do intervalo de confiança para média ($(unit)/ha)", 
            "Limite superior do intervalo de confiança para média ($(unit)/ha)", "Total da população ($(unit))", "Limite inferior do intervalo de confiança para o total ($(unit))", 
            "Limite superior do interva lo de confiança para o total ($(unit))", "Área da população (ha)", "Erro da amostragem relativo (%)", 
            "Erro padrão absoluto ($(unit)/ha)", "Erro padrão da média ($(unit)/ha)", "Variância da média ($(unit)/ha)²", "Fator de correção", 
            "Unidades amostrais possíveis", "Número de unidades amostrais totais", "Número de unidades do inventário florestal", 
            "Número de faixas do inventário florestal",  "Nível de significância (α)"], 
            Valores=round.([mean(Tabela.Média), mean(Tabela.Média)-((quantile(TDist((length(Tabela.n)*first(unique((Tabela.n))))-1), 
            1-α/2))*sqrt(Sx²)), mean(Tabela.Média)+((quantile(TDist((length(Tabela.n)*first(unique((Tabela.n))))-1), 
            1-α/2))*sqrt(Sx²)), ((N*mean(Tabela.Média))/Conversor), (((N*mean(Tabela.Média))-N*((quantile(TDist((length(Tabela.n)*first(unique((Tabela.n))))-1), 
            1-α/2))*sqrt(Sx²)))/Conversor), (((N*mean(Tabela.Média))+N*((quantile(TDist((length(Tabela.n)*first(unique((Tabela.n))))-1), 
            1-α/2))*sqrt(Sx²)))/Conversor), Area, (((quantile(TDist((length(Tabela.n)*first(unique((Tabela.n))))-1), 
            1-α/2))*sqrt(Sx²))/mean(Tabela.Média))*1000, ((quantile(TDist((length(Tabela.n)*first(unique((Tabela.n))))-1), 1-α/2))*sqrt(Sx²)), 
            sqrt(Sx²), Sx², (1-((length(Tabela.n)*first(unique((Tabela.n))))/N)), N, (length(Tabela.n)*first(unique((Tabela.n)))), 
            first(unique((Tabela.n))), length(Tabela.n), α], digits=2))
        ) #Tabela de resultados
    
        return Resultados
    end 
end