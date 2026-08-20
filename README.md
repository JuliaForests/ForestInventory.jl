# ForestInventory

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaForests.github.io/ForestInventory.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaForests.github.io/ForestInventory.jl/dev/)
[![Build Status](https://github.com/JuliaForests/ForestInventory.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/JuliaForests/ForestInventory.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/JuliaForests/ForestInventory.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaForests/ForestInventory.jl)

**ForestInventory.jl** facilitates the processing of forest-inventory sampling
procedures: it estimates a stand-level total (and its confidence interval) — volume,
basal area, tree count, or any other measured quantity — from a set of sample plots or
points, following the sampling designs described in Péllico Netto & Brena's *Inventário
Florestal* (1997). It is the forest-inventory sampling engine behind
[ForestMensuration.jl](https://github.com/JuliaForests/ForestMensuration.jl), and works
standalone too.

## Installation

```julia
using Pkg
Pkg.add("ForestInventory")
```

## Example Usage

Every design takes plain vectors and never assumes a unit for the sampled quantity
itself: give it `u"m^3"`, `u"m^2/ha"`, a plain tree count, whatever it actually is, and
that's exactly what comes back out. `plot_area`/`total_area` take any `Unitful` area
(`ha`, `acre`, `m^2`, ...) and default to hectares when given as bare numbers, matching
the rest of the `JuliaForests` ecosystem.

Estimate a stand total from a simple random sample of plots — the reference design every
other method in the package is compared against. The 70 plot volumes below (1000 m² plots,
a 45 ha stand) come from Nardini's own simple-random-sampling exercise (see [About the
Authors](https://JuliaForests.github.io/ForestInventory.jl/dev/forestinventory/#About-the-Authors)) —
this reproduces its published mean/total (217.97 m³/ha, 9808.71 m³) exactly:

```julia
using ForestInventory, ForestFoundations

v = [
    18.5, 22.8, 15.8, 9.8, 23.5, 20.7, 28.3, 6.9, 28.7, 26.4, 7.6, 10.3, 21.1, 26.3, 39.9,
    12.7, 11.5, 11.0, 18.9, 31.4, 9.6, 9.8, 7.1, 25.9, 25.7, 24.2, 19.0, 25.9, 37.7, 23.3,
    24.6, 32.3, 23.4, 17.1, 8.9, 32.7, 8.1, 29.2, 24.7, 25.6, 13.2, 25.7, 27.0, 20.4, 13.6,
    32.4, 31.3, 34.6, 9.5, 18.7, 28.9, 21.6, 8.3, 23.2, 26.6, 25.6, 8.0, 22.1, 21.9, 23.3,
    9.9, 30.7, 29.2, 28.0, 13.7, 27.6, 37.9, 27.3, 28.0, 30.7,
]u"m^3"

simpleCasualSampling(v, 1000u"m^2", 45u"ha")
```

See the [tutorial](https://JuliaForests.github.io/ForestInventory.jl/dev/tutorial/) for
the same treatment of the other 9 designs.

## Sampling designs

| Function | Design |
|---|---|
| `simpleCasualSampling` | Simple random sampling |
| `stratifiedSampling` | Stratified random sampling |
| `systematicSampling` | Systematic sampling (method of successive differences) |
| `multistartSystematicSampling` | Systematic sampling with multiple random starts |
| `clusterSampling` | One-stage cluster sampling |
| `twoStageSampling` | Two-stage sampling |
| `horizontalPointSampling` | Horizontal point sampling (Bitterlich's angle-count method) |
| `independentOccasionsSampling` | Successive occasions, independent samples |
| `completeReplacementSampling` | Successive occasions, complete replacement (matched plots) |
| `partialReplacementSampling` | Successive occasions, partial replacement |
| `doubleSampling` | Successive occasions, double sampling with regression |

See the [documentation](https://JuliaForests.github.io/ForestInventory.jl/dev/) for a
full worked example of every design.

## About the Authors

This package continues a lineage of academic work on forest-inventory software written in
Julia.

It began as a desktop application built for a Master's dissertation:

> Nardini, C. (2023). *Desenvolvimento de Software para Uso em Inventário Florestal
> Baseado em Linguagem Julia*. Dissertação (Mestrado em Agronomia: Agricultura e Ambiente)
> — Universidade Federal de Santa Maria, Campus Frederico Westphalen.
> [PDF](https://repositorio.ufsm.br/bitstream/handle/1/28263/DIS_PPGAAA_2023_NARDINI_CLAITON.pdf?sequence=1&isAllowed=y#page=72)

That application — [ForestInventory.jl](https://github.com/renilsonlisboa/ForestInventory.jl)
by Renilson Lisboa Júnior, itself built on Péllico Netto & Brena's *Inventário Florestal*
(1997) — was later taken further as an undergraduate final course project (TCC), which
turned it into a web platform:

> Lisboa Júnior, R. (2024). *Plataforma Web para Processamento de Inventário Florestal em
> Julia Language*. Trabalho de Conclusão de Curso (Bacharel em Engenharia Florestal) —
> Universidade Federal do Paraná, Curitiba.
> [PDF](https://creaweb3.crea-pr.org.br/webcrea/biblioteca_virtual/arquivos_BV/2025/36_TCC-Renilson.pdf)

This package, **ForestInventory.jl**, rewrites and reinterprets that same set of 11
sampling designs by Marcos Daniel da Silva, a Forest Engineer, integrating them into the
`JuliaForests` ecosystem as part of his own Master's dissertation, under the guidance of
Prof. Dr. Emanuel Araújo da Silva (advisor) and Prof. Rafaelo Balbinot (co-advisor), in
the Programa de Pós-Graduação em Agronomia, Agricultura e Ambiente (PPGAAA) at the
Universidade Federal de Santa Maria, Campus Frederico Westphalen (UFSM-FW):

> Silva, M. D. da. *JuliaForests: da teoria à análise computacional de Florestas*.
> Dissertação (Mestrado em Agronomia: Agricultura e Ambiente) — Universidade Federal de
> Santa Maria, Campus Frederico Westphalen. Orientador: Prof. Dr. Emanuel Araújo da Silva;
> coorientador: Prof. Rafaelo Balbinot.

Contact: [LinkedIn](https://www.linkedin.com/in/marcosdanieldasilva/?locale=en_US) ·
[marcosdasilva@5a.tec.br](mailto:marcosdasilva@5a.tec.br).

## License

MIT
