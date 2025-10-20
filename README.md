README bilíngue/Bilingual README

******************************************************************************
# Momentum Transport Analysis in the TCABR Tokamak

This repository contains the MATLAB source code for the analysis of toroidal momentum transport in the TCABR tokamak. This work is a modular and robust refactoring of the original code used in the doctoral thesis, "A Study of Momentum Transport from the Toroidal Rotation Velocity of Plasma and its Interaction with Neutral Particles in the TCABR Tokamak".

The primary objective of this codebase is to calculate the effective momentum diffusivity ($\chi_{\varphi}^{eff}$) based on a method that utilises a Fourier-Bessel series expansion of the experimental velocity profile and the ion-neutral collision frequency. The results are then compared with theoretical transport models from the literature.

## Key Features

* **Statistical Profile Analysis:** Performs a robust analysis of experimental toroidal velocity ($V_\varphi$) and ion temperature ($T_i$) profiles using a bootstrap Monte Carlo method to quantify uncertainties.
* **Fourier-Bessel Series Expansion:** Fits the mean velocity profile to a Bessel function series to determine the system's eigenvalues ($\lambda_j$), which are fundamental to the diffusivity calculation.
* **Derived Profile Calculation:** Computes a comprehensive set of derived physical profiles, including the safety factor (q), magnetic shear (s), collision times and frequencies, and collisionality ($\nu_*$).
* **Effective Diffusivity Calculation:** Implements the core method from the thesis to calculate the $\chi_{\varphi}^{eff}$ profile from the ion-neutral collision frequency and the derived eigenvalues.
* **Theoretical Models:** Calculates transport profiles based on theoretical models from the literature, such as the Helander model for velocity and the Solomon model for diffusivity and pinch velocity, for comparison purposes.
* **Caching System:** Includes a system to save and load intermediate results, avoiding the need to recalculate everything upon each execution.

## Project Structure

The code is organised into a modular structure to ensure clarity and maintainability.

/tcabr_momentum_codes/ 
|-- main_analysis.m # Main script that orchestrates the entire analysis. 
| |-- /src/ # Contains all calculation and analysis functions. 
| |-- analyze_velocity_profile.m 
| |-- analyze_temperature_profile.m 
| |-- compute_effective_diffusivity.m 
| |-- ... (and others) 
| |-- /plotting/ # Contains dedicated functions for generating plots. 
| |-- run_verification_plots.m 
| |-- plot_diffusivity_comparison.m 
| |-- /utils/ # Auxiliary functions (e.g., findBesselZeros.m). 
| |-- /data/ # Input data files (.txt). 
| |-- /results/ # Directory for outputs (plots, .mat files). 
| |-- README.md # This fil

## Requirements

* **MATLAB** (version R2021a or newer recommended)
* **Parallel Computing Toolbox™:** Required for the accelerated execution of `parfor` loops.
* **Optimisation Toolbox™:** Required for the `lsqnonlin` function used in the temperature analysis.
* **Curve Fitting Toolbox™:** Required for the `fit` function used in the velocity analysis.

## How to Run

1.  Clone this repository to your local machine.
2.  Open MATLAB.
3.  Navigate to the project's root folder (`/tcabr_momentum_codes/`).
4.  Execute the main script from the Command Window:
    ```matlab
    >> main_analysis.m
    ```

### Recalculation Control

The `main_analysis.m` script features a switch at the top named `force_recalculation`.
* If `force_recalculation = false;` (default), the script will attempt to load the results of a previous run from the `/results/full_analysis_results.mat` file, skipping all computationally intensive calculations.
* If `force_recalculation = true;`, the script will run all calculations from scratch and save the new results at the end. Change this to `true` whenever you make modifications to any function in the `/src/` folder.

## Methodology Overview

The core method implemented for the calculation of the effective diffusivity follows the derivation presented in the thesis, culminating in **Equation (3.13)**:
$$\chi_{\varphi}^{eff} = \sum_{j}\frac{1}{\lambda_{j}^{2}} \left( \frac{3}{4\epsilon} - 1 \right) \nu_{iH_{0}}$$

* The eigenvalues **$\lambda_j$** are determined from the zeros of a Fourier-Bessel series expansion of the experimental toroidal velocity profile ($V_{\varphi}$).
* The ion-neutral collision frequency **$\nu_{iH_{0}}$** is calculated as $\nu_{iH_{0}} = \langle\sigma v\rangle_{cx} n_{H_{0}}(r)$, where the charge-exchange rate coefficient $\langle\sigma v\rangle_{cx}$ depends on the ion temperature profile ($T_i$) and $n_{H_{0}}$ is the neutral density profile.

## Licence

This project is licensed under the **MIT Licence**. See the `LICENSE` file for further details.

## How to Cite

If you use this code in your research, please cite the original thesis:

> Novaes, D. O. (2025). *A Study of Momentum Transport from the Toroidal Rotation Velocity of Plasma and its Interaction with Neutral Particles in the TCABR Tokamak*. Doctoral Thesis, Postgraduate Programme in Physics, Federal University of Rio Grande do Sul, Porto Alegre, Brazil.

## Author and Contact

* **Douglas Oliveira Novaes**
* Email: `dougnovaes@alumni.usp.br`


******************************************************************************
# Análise de Transporte de Momento no Tokamak TCABR

Este repositório contém o código-fonte em MATLAB para a análise do transporte de momento toroidal no tokamak TCABR. O trabalho é uma refatoração modular e robusta do código original utilizado na tese de doutorado "Estudo do Transporte de Momento a Partir da Velocidade de Rotação Toroidal de Plasma em sua Interação com Partículas Neutras no Tokamak TCABR".

O objetivo principal deste código é calcular a difusividade de momento efetiva ($\chi_{\varphi}^{eff}$) com base em um método que utiliza a expansão em séries de Fourier-Bessel do perfil de velocidade experimental e a frequência de colisão íon-neutro. Os resultados são então comparados com modelos teóricos de transporte da literatura.

## Principais Funcionalidades

* **Análise Estatística de Perfis:** Realiza uma análise robusta dos perfis experimentais de velocidade toroidal ($V_\varphi$) e temperatura iônica ($T_i$) usando um método de bootstrap Monte Carlo para quantificar incertezas.
* **Expansão em Séries de Fourier-Bessel:** Ajusta o perfil de velocidade médio a uma série de funções de Bessel para determinar os autovalores ($\lambda_j$) do sistema, que são fundamentais para o cálculo da difusividade.
* **Cálculo de Perfis Derivados:** Calcula um conjunto completo de perfis físicos derivados, incluindo fator de segurança (q), cisalhamento magnético (s), tempos e frequências de colisão, e colisionalidade ($\nu_*$).
* **Cálculo da Difusividade Efetiva:** Implementa o método central da tese para calcular o perfil de $\chi_{\varphi}^{eff}$ a partir da frequência de colisão íon-neutro e dos autovalores.
* **Modelos Teóricos:** Calcula perfis de transporte baseados em modelos teóricos da literatura, como o modelo de Helander para a velocidade e o modelo de Solomon para a difusividade e velocidade de pinch, para fins de comparação.
* **Sistema de Cache:** Inclui um sistema para salvar e carregar resultados intermediários, evitando a necessidade de recalcular tudo a cada execução.

## Estrutura do Projeto

O código é organizado em uma estrutura modular para garantir clareza e manutenibilidade.

/tcabr_momentum_codes/ 
|-- main_analysis.m # Script principal que orquestra toda a análise. 
| |-- /src/ # Contém todas as funções de cálculo e análise. 
| |-- analyze_velocity_profile.m 
| |-- analyze_temperature_profile.m 
| |-- compute_effective_diffusivity.m 
| |-- ... (e outras) 
| |-- /plotting/ # Contém funções dedicadas à geração de gráficos. 
| |-- run_verification_plots.m 
| |-- plot_diffusivity_comparison.m 
| |-- /utils/ # Funções auxiliares (ex: findBesselZeros.m). 
| |-- /data/ # Arquivos de dados de entrada (.txt). 
| |-- /results/ # Diretório para saídas (gráficos, arquivos .mat). 
| |-- README.md # Este arquivo.

## Requisitos

* **MATLAB** (versão R2021a ou mais recente recomendada)
* **Parallel Computing Toolbox™:** Necessário para a execução acelerada dos laços `parfor`.
* **Optimization Toolbox™:** Necessário para a função `lsqnonlin` usada na análise de temperatura.
* **Curve Fitting Toolbox™:** Necessário para a função `fit` usada na análise de velocidade.

## Como Executar

1.  Clone este repositório para a sua máquina local.
2.  Abra o MATLAB.
3.  Navegue até a pasta raiz do projeto (`/tcabr_momentum_codes/`).
4.  Execute o script principal a partir da Janela de Comando:
    ```matlab
    >> main_analysis.m
    ```

### Controle de Recálculo

O `main_analysis.m` possui um "interruptor" na parte superior chamado `force_recalculation`.
* Se `force_recalculation = false;` (padrão), o script tentará carregar os resultados de uma execução anterior do arquivo `/results/full_analysis_results.mat`, pulando todos os cálculos pesados.
* Se `force_recalculation = true;`, o script rodará todos os cálculos desde o início e salvará os novos resultados no final. Mude para `true` sempre que fizer alterações em alguma função na pasta `/src/`.

## Visão Geral da Metodologia

O método central implementado para o cálculo da difusividade efetiva segue a derivação apresentada na tese, culminando na **Equação (3.13)**:
$$\chi_{\varphi}^{eff} = \sum_{j}\frac{1}{\lambda_{j}^{2}} \left( \frac{3}{4\epsilon} - 1 \right) \nu_{iH_{0}}$$

* Os autovalores **$\lambda_j$** são determinados a partir dos zeros de uma expansão em série de Fourier-Bessel do perfil de velocidade toroidal experimental ($V_{\varphi}$).
* A frequência de colisão íon-neutro **$\nu_{iH_{0}}$** é calculada como $\nu_{iH_{0}} = \langle\sigma v\rangle_{cx} n_{H_{0}}(r)$, onde o coeficiente de taxa de troca de carga $\langle\sigma v\rangle_{cx}$ depende do perfil de temperatura iônica ($T_i$) e $n_{H_{0}}$ é o perfil de densidade de neutros.

## Licença

Este projeto é licenciado sob a **Licença MIT**. Veja o arquivo `LICENSE` para mais detalhes.

## Como Citar

Se você utilizar este código em seu trabalho de pesquisa, por favor, cite a tese original:

> Novaes, D. O. (2025). *Estudo do Transporte de Momento a Partir da Velocidade de Rotação Toroidal de Plasma em sua Interação com Partículas Neutras no Tokamak TCABR*. Tese de Doutorado, Programa de Pós-Graduação em Física, Universidade Federal do Rio Grande do Sul, Porto Alegre, Brasil.

## Autor e Contato

* **Douglas Oliveira Novaes**
* E-mail: `dougnovaes@alumni.usp.br`
***********************************************************************

