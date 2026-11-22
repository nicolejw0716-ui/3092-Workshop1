Enzyme Inhibitor Lab

Fit your kinetics data. See what your inhibitor is doing.

The idea
Enzyme Inhibitor Lab is a browser app. You upload enzyme kinetics measurements
(substrate concentration and reaction velocity) for a control and one or more
inhibitor conditions. For each condition, the app:

1. fits the Michaelis–Menten equation, v = Vmax·[S] / (Km + [S]), with R's `nls()`
2. reports Vmax and Km with standard errors
3. plots the raw measurements with the fitted curves and dashed Vmax lines
4. compares each condition with the control and suggests the type of inhibition
   (competitive, non-competitive, mixed)

It runs on webR (R compiled to WebAssembly), so all calculations happen in
the browser and no data leaves your device.

## Link to biochemistry research
Characterising inhibitors is a standard step in drug discovery and enzymology.
The type of inhibition says something about the mechanism:

| Type | Vmax | Km | Meaning |
|---|---|---|---|
| Competitive | unchanged | increases | inhibitor competes with substrate for the active site |
| Non-competitive | decreases | unchanged | inhibitor binds elsewhere, equally to E and ES |
| Mixed | decreases | changes | inhibitor binds E and ES with different affinities |

The app fits the nonlinear equation directly instead of using a Lineweaver–Burk
plot. Taking reciprocals distorts the experimental error, especially at low
substrate concentrations, so a direct nonlinear fit gives more reliable estimates.

## How to use
1. Prepare a CSV with the columns `condition`, `substrate_mM`, `velocity`
   (one row per measurement; substrate in mM, velocity in µmol/min).
   The first condition is treated as the control.
2. Drop the file onto the upload box or click browse.
3. Click Run analysis.

## Screenshots

Start screen
![alt text](<螢幕截圖 2026-11-22 上午11.17.20.png>)

Results
![alt text](<螢幕截圖 2026-11-22 上午11.17.56.png>)