# ForestInventory Package:

Documentation for [ForestInventory](https://github.com/JuliaForests/ForestInventory.jl).

**ForestInventory.jl** facilitates the processing of forest-inventory sampling procedures in Julia: it estimates a stand-level total (and its confidence interval) — volume, basal area, tree count, or any other measured quantity — from a set of sample plots or points, following the sampling designs described in Péllico Netto & Brena's *Inventário Florestal* (1997). It is the forest-inventory sampling engine behind [ForestMensuration.jl](https://github.com/JuliaForests/ForestMensuration.jl), and works standalone too. Key features include:

- **11 Classic Sampling Designs**: Simple, stratified, systematic (with an optional multiple-random-start variant), one-stage cluster, two-stage, horizontal point/Bitterlich angle-count, and the four sampling-on-successive-occasions designs used in continuous forest inventory — each its own function, generic to any number of strata/clusters/plots.
- **Plain-Vector Interface**: Every design takes plain vectors — no `DataFrame`/column-name indirection — so it drops straight into any existing data-handling pipeline.
- **Unit-Aware, But Not Unit-Assuming**: Built on `Unitful`, like the rest of the `JuliaForests` ecosystem, but the sampled quantity itself is never assumed to be volume — a `Unitful` quantity is preserved exactly as given, and a plain number stays unitless all the way through.
- **Structured Reports**: Designs that produce several related tables (one per stratum/cluster, or one per inventory occasion) return a [`SamplingReport`](@ref), whose tables are reachable through accessor functions or directly as properties.

Whether you are running a single simple-random-sample estimate or a full continuous forest inventory across successive occasions, **ForestInventory.jl** simplifies the statistics so you can focus on the fieldwork and the resulting decisions.

```@meta
CurrentModule = ForestInventory
```

## Installation

The ForestInventory package is available through the Julia package system and can be added by running Pkg.add("ForestInventory") or by directly downloading it from the GitHub page:

```julia-repl
pkg> add https://github.com/JuliaForests/ForestInventory.jl
```

## About the Authors

This package continues a lineage of academic work on forest-inventory software written in Julia.

It began as a desktop application built for a Master's dissertation:

> Nardini, C. (2023). *Desenvolvimento de Software para Uso em Inventário Florestal Baseado em Linguagem Julia*. Dissertação (Mestrado em Agronomia: Agricultura e Ambiente) — Universidade Federal de Santa Maria, Campus Frederico Westphalen. Orientador: Braulio Otomar Caron; coorientador: Cristiano Bertolini. [PDF](https://repositorio.ufsm.br/bitstream/handle/1/28263/DIS_PPGAAA_2023_NARDINI_CLAITON.pdf?sequence=1&isAllowed=y#page=72)

That application — [ForestInventory.jl](https://github.com/renilsonlisboa/ForestInventory.jl) by Renilson Lisboa Júnior, itself built on Péllico Netto & Brena's *Inventário Florestal* (1997) — was later taken further as an undergraduate final course project (TCC), which turned it into a web platform:

> Lisboa Júnior, R. (2024). *Plataforma Web para Processamento de Inventário Florestal em Julia Language*. Trabalho de Conclusão de Curso (Bacharel em Engenharia Florestal) — Universidade Federal do Paraná, Curitiba. Orientador: Alexandre Behling; coorientadores: Richardson Ribeiro, Claiton Nardini. [PDF](https://creaweb3.crea-pr.org.br/webcrea/biblioteca_virtual/arquivos_BV/2025/36_TCC-Renilson.pdf)

This package, **ForestInventory.jl**, is a from-scratch, Julia-idiomatic reimplementation of that same set of 11 sampling designs — proper package structure, `Unitful` support throughout, a full test suite, and generated documentation — rather than a fork or a port of the original GUI application. I am Marcos Daniel da Silva, a Forest Engineer, and this reimplementation — together with the rest of the `JuliaForests` ecosystem ([ForestMensuration.jl](https://github.com/JuliaForests/ForestMensuration.jl), [ForestModeling.jl](https://github.com/JuliaForests/ForestModeling.jl), [ForestFoundations.jl](https://github.com/JuliaForests/ForestFoundations.jl)) — is developed as part of my own Master's dissertation, under the guidance of Prof. Dr. Emanuel Araújo da Silva (advisor) and Prof. Rafaelo Balbinot (co-advisor), in the Programa de Pós-Graduação em Agronomia, Agricultura e Ambiente (PPGAAA) at the Universidade Federal de Santa Maria, Campus Frederico Westphalen (UFSM-FW):

> Silva, M. D. da. *JuliaForests: da teoria à análise computacional de Florestas*. Dissertação (Mestrado em Agronomia: Agricultura e Ambiente) — Universidade Federal de Santa Maria, Campus Frederico Westphalen. Orientador: Prof. Dr. Emanuel Araújo da Silva; coorientador: Prof. Rafaelo Balbinot.

For more details, please feel free to connect with me on [LinkedIn](https://www.linkedin.com/in/marcosdanieldasilva/?locale=en_US) or by e-mail: [marcosdasilva@5a.tec.br](mailto:marcosdasilva@5a.tec.br).

This package is supported and encouraged by [5A Inteligência e Engenharia](https://5a.tec.br/).
