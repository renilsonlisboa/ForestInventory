# ForestInventory

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://renilsonlisboa.github.io/ForestInventory.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://renilsonlisboa.github.io/ForestInventory.jl/dev/)
[![Build Status](https://github.com/renilsonlisboa/ForestInventory.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/renilsonlisboa/ForestInventory.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/renilsonlisboa/ForestInventory.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/renilsonlisboa/ForestInventory.jl)

## ForestInventory

<p align="justify"> O pacote ForestInventory.jl é um aplicativo desenvolvido em Linguagem Julia. Ele foi projetado para facilitar os processamentos dos processos de amostragem do inventário florestal. Nele é possível analisar 10 processos de amostragem do inventário florestal, que estão baseados no Livro de Inventário Florestal, 1997 de Sylvio Péllico Netto e Doádi Antônio Brena. Essa ferramenta tem o intuito de contribuir e facilitar o processamento de dados por pesquisadores, engenheiros florestais, estatísticos e pessoas que trabalhem com dados relacionados ao inventário florestal. Além disso, o seu uso é de livre acesso, possibilitando seu uso no ensino superior.</p>

## Métodos de Amostragem

- [Área Fixa](https://repositorio.ufsm.br/bitstream/handle/1/28263/DIS_PPGAAA_2023_NARDINI_CLAITON.pdf?sequence=1&isAllowed=y#page=60)
- [Bitterlich](https://repositorio.ufsm.br/bitstream/handle/1/28263/DIS_PPGAAA_2023_NARDINI_CLAITON.pdf?sequence=1&isAllowed=y#page=63)
- [Prodan](https://repositorio.ufsm.br/bitstream/handle/1/28263/DIS_PPGAAA_2023_NARDINI_CLAITON.pdf?sequence=1&isAllowed=y#page=66)
- [Strand](https://repositorio.ufsm.br/bitstream/handle/1/28263/DIS_PPGAAA_2023_NARDINI_CLAITON.pdf?sequence=1&isAllowed=y#page=60)
- [3-P](https://repositorio.ufsm.br/bitstream/handle/1/28263/DIS_PPGAAA_2023_NARDINI_CLAITON.pdf?sequence=1&isAllowed=y#page=6)


## Processos de Amostragem

- [Amostragem Aleatória Simples](https://repositorio.ufsm.br/bitstream/handle/1/28263/DIS_PPGAAA_2023_NARDINI_CLAITON.pdf?sequence=1&isAllowed=y#page=60)
- [Amostragem Estratificada](https://repositorio.ufsm.br/bitstream/handle/1/28263/DIS_PPGAAA_2023_NARDINI_CLAITON.pdf?sequence=1&isAllowed=y#page=63)
- [Amostragem Sistemática](https://repositorio.ufsm.br/bitstream/handle/1/28263/DIS_PPGAAA_2023_NARDINI_CLAITON.pdf?sequence=1&isAllowed=y#page=66)
- [Amostragem em Dois Estágios](https://repositorio.ufsm.br/bitstream/handle/1/28263/DIS_PPGAAA_2023_NARDINI_CLAITON.pdf?sequence=1&isAllowed=y#page=69)
- [Amostragem em Conglemerados](https://repositorio.ufsm.br/bitstream/handle/1/28263/DIS_PPGAAA_2023_NARDINI_CLAITON.pdf?sequence=1&isAllowed=y#page=72)
- [Amostragem Sistemática com Múltiplos Inícios Alatórios](https://repositorio.ufsm.br/bitstream/handle/1/28263/DIS_PPGAAA_2023_NARDINI_CLAITON.pdf?sequence=1&isAllowed=y#page=75)
- [Amostragem Independente](https://repositorio.ufsm.br/bitstream/handle/1/28263/DIS_PPGAAA_2023_NARDINI_CLAITON.pdf?sequence=1&isAllowed=y#page=78)
- [Amostragem com Repetição Total](https://repositorio.ufsm.br/bitstream/handle/1/28263/DIS_PPGAAA_2023_NARDINI_CLAITON.pdf?sequence=1&isAllowed=y#page=81)
- [Amostragem com Repetição Dupla](https://repositorio.ufsm.br/bitstream/handle/1/28263/DIS_PPGAAA_2023_NARDINI_CLAITON.pdf?sequence=1&isAllowed=y#page=84)
- [Amostragem com Repetição Parcial](https://repositorio.ufsm.br/bitstream/handle/1/28263/DIS_PPGAAA_2023_NARDINI_CLAITON.pdf?sequence=1&isAllowed=y#page=87)

## Instalação do Pacote

<p align="justify"> A instalação do pacote ForesteInventory.jl e suas dependências pode ser feita por meio do gestor de pacotes de Julia. Para isto, basta o usuário abrir o Julia REPL e executar os comandos: </p>

```julia
using Pkg

Pkg.add("https://github.com/renilsonlisboa/ForestInventory.jl")
```

## Usando ForestInventory.jl

<p align="justify"> Após finalizada a instalação do pacote por meio dos comandos a cima, o usuário podera utilizar o pacote ForestInventory.jl por meio dos seguintes comandos: </p>

```julia
using ForestInventory

Inventory()
```