# Figures

One folder per manuscript figure, named `fig<N>_<short-name>/` (main figures) or
`figS<N>_<short-name>/` (supplementary). Each figure folder reads its inputs from
[`../results/tables`](../results/tables) (written by [`../analysis`](../analysis)) and/or the
processed objects in `../data/processed` (written by [`../preprocessing`](../preprocessing)) —
figure scripts should not re-run preprocessing or analysis themselves.

Suggested contents per figure folder:
- one script (`.R` or `.py`) per figure/panel that produces the final plot(s)
- output plots written to `../results/figures/`

This folder is currently a scaffold — no figure scripts have been added yet.
