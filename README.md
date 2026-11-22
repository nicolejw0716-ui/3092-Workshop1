# Enzyme Inhibitor Lab

A browser-based enzyme kinetics workbench. Upload a CSV with `condition`,
`substrate_mM`, and `velocity` columns to fit a separate Michaelis–Menten
`nls()` model for each condition using webR and base R. The example dataset is
loaded automatically from `data/enzyme_kinetics.csv`.

## Run locally

```sh
npm install
npm run dev
```

The analysis runs in WebAssembly in your browser; uploaded measurements are
not sent to a server. Rows missing a condition, containing non-numeric or
non-finite measurements, or with a negative substrate concentration are
reported as dropped and excluded from the fits.

Results include a base R plot and Vmax/Km estimates with standard errors
rounded to two significant figures. Interpretation compares each condition
with `control` (or the first condition), treating parameter differences within
one combined standard error as similar; changes in both parameters are labeled
consistent with mixed inhibition.