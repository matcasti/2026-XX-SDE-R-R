# SDE-Inverse Gaussian HRV Framework

Reproducible research compendium for:

> **A Continuous-Time Stochastic State-Space Framework for Heart Rate Variability: The SDE-Inverse Gaussian Model**
>
> Matías Castillo-Aguilar — Escuela de Medicina, Universidad de Magallanes

---

## Repository structure

```
.
├── index.qmd              # Main manuscript (Quarto)
├── _quarto.yml            # Quarto project config
├── styles.css             # Custom CSS
├── references.bib         # BibTeX references
├── renv_setup.R           # One-time environment setup (run first)
├── renv.lock              # Pinned package versions (committed)
│
└── R/
    ├── model_functions.R  # Canonical generative model (sourced by index.qmd)
    ├── mle.R              # Analytic MLE and FIM functions
    └── identifiability.R  # Recovery study and profile likelihood
```

---

## Reproducing the paper

### Requirements

- R ≥ 4.3.0
- Quarto ≥ 1.4
- No external system dependencies beyond a standard R installation

### Step 1 — Set up the R environment (run once)

```bash
Rscript renv_setup.R
```

This installs all required packages at pinned versions and writes `renv.lock`.

### Step 2 — Render the manuscript

```bash
quarto render index.qmd
```

This single command generates all figures, tables, and numerical results in
the paper. No manual steps, no external data files, no pre-computed results.
The document is fully self-contained.

### Step 3 — Verify reproducibility

All random seeds are set explicitly:

| Seed | Location | Purpose |
|------|----------|---------|
| `2026` | `index.qmd` setup chunk | Main simulation (`RES`) |
| `101`  | `index.qmd` latent-state figure | Baseline OU visualization |
| `123`  | `index.qmd` log-link figure | Log vs linear link comparison |
| `321`  | `index.qmd` OU bias figure | Euler-Maruyama comparison |
| `2026` | `identifiability.R` (`RECOVERY_MASTER_SEED`) | Recovery study seed stream |

Rendering on a different platform may produce bit-identical results only
within the same R version and OS due to `.Random.seed` portability.
Numerical results should be stable to the precision reported in the paper.

---

## Computational notes

| Component | Approx. runtime (single core) | Notes |
|-----------|-------------------------------|-------|
| Main simulation (720 s, dt=0.005) | ~15 s | Cached after first render |
| Recovery study (N=200, 300 s each) | ~5 min | Parallelized on Unix; cached |
| Profile likelihoods (6 parameters) | ~30 s | Cached |
| All figures and tables | ~2 min | — |

The recovery study uses `parallel::mclapply` on Unix/macOS (falls back to
`lapply` on Windows). The Quarto `cache: true` option on expensive chunks
avoids re-running on subsequent renders if inputs are unchanged.

---

## Archival

The rendered HTML output, code, and `renv.lock` are archived at:

- **OSF / Zenodo**: [DOI to be assigned upon submission]

To regenerate the archived output from scratch:

```bash
Rscript renv_setup.R   # restore pinned environment
quarto render index.qmd
```

---

## License

Code: MIT. Manuscript text: CC BY 4.0.
