# Physics-Based Estimation of Toroidal Momentum Transport in Tokamaks

This repository implements a physics-based method to estimate the **effective toroidal momentum diffusivity** ($χ_ϕ^\text{eff}$) in tokamak plasmas, with explicit treatment of **ion–neutral interactions**, applied to experimental rotation profiles from the **TCABR tokamak**.

The method relies on a Fourier–Bessel expansion of the measured toroidal velocity profile and on the ion–neutral collision frequency, enabling a direct estimation of $χ_ϕ^\text{eff}$ in plasmas with **no external momentum input**. The results are strictly validated against **collisionality regimes** (ν\*) and compared with established neoclassical and semi-empirical transport models.

This code accompanies and reproduces the analysis presented in:

> **Novaes, D. O.**, Severo, J. H. F., Rizzato, F. B. et al.  
> *Estimation of Effective Momentum Diffusivity and Its Correlation with Neutral Particle Density Based on Toroidal Rotation Profiles in the TCABR Tokamak*.  
> **Brazilian Journal of Physics**, 55, 43 (2025).  
> DOI: [10.1007/s13538-024-01681-x](https://doi.org/10.1007/s13538-024-01681-x)

---

## Scientific Scope

### Primary Output
- Radial profile of the effective toroidal momentum diffusivity $χ_ϕ^\text{eff}$.
- Quantitative comparison with neoclassical (Helander) and semi-empirical transport models.
- Collisionality-based regime validation (Pfirsch–Schlüter).

### What this code is
- A reproducible analysis pipeline linking experimental rotation data to momentum transport coefficients.
- A tool for studying intrinsic rotation and ion–neutral momentum transport in medium-size tokamaks.

### What this code is not
- A general-purpose transport solver.
- A predictive integrated modelling framework.

---

## Method Summary

The effective momentum diffusivity is computed following the derivation detailed in the reference article. The core estimator is:

$\chi_{\varphi}^\text{eff} = \sum_{j}\frac{1}{\lambda_{j}^{2}}
\left( \frac{3}{4\epsilon} - 1 \right)\nu_\text{iH₀}$

where:
- the eigenvalues λⱼ, from a Helmholtz-like equation, are obtained from the zeros of a Fourier–Bessel expansion fitted to the experimental toroidal velocity profile,
- the ion–neutral collision frequency $ν_\text{iH₀}$ depends on the ion temperature and neutral density profiles,
- experimental uncertainties are propagated using a Monte Carlo approach.

---

## Code Structure

- `main_analysis.m` — main analysis pipeline.
- `src/` — core physics and analysis routines.
- `plotting/` — plotting utilities.
- `data/` — experimental input profiles.
- `results/` — cached outputs.

---

## Requirements

- MATLAB R2021a or newer
- Parallel Computing Toolbox
- Optimization Toolbox
- Curve Fitting Toolbox

---

## Licence

MIT Licence.

---

## Author

**Douglas Oliveira Novaes**  
dougnovaes@alumni.usp.br
