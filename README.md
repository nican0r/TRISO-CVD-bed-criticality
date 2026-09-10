# TRISO-CVD-bed-criticality

OpenMC criticality safety (NCS) screening model of a laboratory-scale fluidized-bed chemical vapor deposition furnace used to coat TRISO fuel kernels. The model evaluates whether HALEU UCO TRISO particles accumulating inside a graphite retort during the CVD coating process can approach criticality under normal and off-normal (flooding) conditions. All results are unvalidated screening calculations until benchmarked against ICSBEP evaluations.

## Steps

| Step | Summary | Details |
|------|---------|---------|
| Preamble | Requirements, working agreement, physical description, ambiguity resolution | [docs/steps/preamble.md](docs/steps/preamble.md) |
| 0 — Scaffold | Repo layout, `params.yaml`, `load_params()`, `check_env()` | [docs/steps/step-0-scaffold.md](docs/steps/step-0-scaffold.md) |
| 1 — Materials | All OpenMC material definitions (kernel, coatings, graphite, gas, water, air) | [docs/steps/step-1-materials.md](docs/steps/step-1-materials.md) |
| 2 — TRISO packed bed | Single-particle universe factory, random bed packing with disk cache, lattice optimisation, bed statistics, analytical stage-progression table | [docs/steps/step-2-triso-packed-bed.md](docs/steps/step-2-triso-packed-bed.md) |
| 3 — Bed and cone geometry | Hybrid exact-cone/staircase bed geometry, two bed states (fluidized/collapsed), mass-conserving packing, n_slabs discretization convergence study | [docs/steps/step-3-bed-cone-geometry.md](docs/steps/step-3-bed-cone-geometry.md) |
| 4 — Retort wall, heater, injector | Full furnace shell: graphite retort and cone walls, vacuum gap, heater annulus, water-cooled injector with annular coolant channel, outer vacuum boundary | [docs/steps/step-4-retort-heater-injector.md](docs/steps/step-4-retort-heater-injector.md) |
