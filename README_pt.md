# Estimativa Física do Transporte de Momento Toroidal em Tokamaks

Este repositório implementa um método baseado em princípios físicos para a estimativa da **difusividade efetiva de momento toroidal** ($\chi_{\varphi}^\text{eff}$) em plasmas de tokamak, com tratamento explícito das **interações íon–neutro**, aplicado a perfis experimentais de rotação do **tokamak TCABR**.

O método baseia-se na expansão em séries de Fourier–Bessel do perfil experimental de velocidade toroidal e na frequência de colisão íon–neutro, permitindo a estimativa direta de $\chi_{\varphi}^\text{eff}$ em plasmas **sem injeção externa de momento**. Os resultados são rigorosamente validados em função dos **regimes de colisionalidade** (ν\*) e comparados com modelos neoclássicos e semiempíricos de transporte.

Este código acompanha e reproduz a análise apresentada em:

> **Novaes, D. O.**, Severo, J. H. F., Rizzato, F. B. et al.  
> *Estimation of Effective Momentum Diffusivity and Its Correlation with Neutral Particle Density Based on Toroidal Rotation Profiles in the TCABR Tokamak*.  
> **Brazilian Journal of Physics**, 55, 43 (2025).  
> DOI: 10.1007/s13538-024-01681-x

---

## Escopo Científico

### Resultados Principais
- Perfil radial da difusividade efetiva de momento toroidal $\chi_{\varphi}^\text{eff}$.
- Comparação quantitativa com modelos neoclássicos (Helander) e semiempíricos de transporte.
- Validação do regime físico com base na colisionalidade (Pfirsch–Schlüter).

### O que este código é
- Um pipeline reprodutível que conecta dados experimentais de rotação a coeficientes de transporte de momento.
- Uma ferramenta para o estudo da rotação intrínseca e do transporte de momento via interações íon–neutro em tokamaks de médio porte.

### O que este código não é
- Um solucionador geral de transporte.
- Uma ferramenta de modelagem integrada preditiva.

---

## Resumo do Método

A difusividade efetiva de momento é calculada a partir da derivação apresentada no artigo de referência. O estimador central é dado por:

$$
\chi_{\varphi}^\text{eff} = \sum_{j}\frac{1}{\lambda_{j}^{2}}
\left( \frac{3}{4\epsilon} - 1 \right)\nu_\text{iH₀}
$$

onde:
- os autovalores λⱼ são obtidos a partir dos zeros da expansão em série de Fourier–Bessel ajustada ao perfil experimental de velocidade toroidal,
- a frequência de colisão íon–neutro $\nu_\text{iH₀}$ depende dos perfis de temperatura iônica e densidade de neutros,
- as incertezas experimentais são propagadas por meio de um método de bootstrap Monte Carlo.

---

## Estrutura do Código

- `main_analysis.m` — pipeline principal da análise.
- `src/` — rotinas centrais de física e análise.
- `plotting/` — utilitários de geração de gráficos.
- `data/` — perfis experimentais de entrada.
- `results/` — resultados em cache.

---

## Requisitos

- MATLAB R2021a ou mais recente
- Parallel Computing Toolbox
- Optimization Toolbox
- Curve Fitting Toolbox

---

## Licença

Licença MIT.

---

## Autor

**Douglas Oliveira Novaes**  
dougnovaes@alumni.usp.br
