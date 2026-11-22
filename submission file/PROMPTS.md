GOAL: Build a Vite web app called "Honest Lab Graph Maker" that runs R in the browser using webR.
DATA: The user uploads a CSV with columns: group, replicate, and one numeric measurement column. Load data/enzyme_activity.csv as the default example. Each row is one biological replicate.
METHOD: For each group calculate mean, SD and n. If there are 2 groups, run t.test(); if 3 or more groups, run aov() followed by TukeyHSD(). Do not run multiple separate t-tests.
OUTPUT: (1) A plot drawn with base R graphics showing every raw data point, plus the mean ± SD for each group. The y-axis must start at 0, and the axis labels must include units. (2) A summary table of mean, SD and n per group. (3) The statistical test results.
CHECK: Show the number of rows read per group, and display the R code that was run under the results.
Then start the app with npm run dev.
## Correction prompt 1
Fix these issues:
1. Use the same colour for each condition in the table and the plot.
2. Make the y-axis start at 0 and extend above the highest fitted Vmax. Draw the fitted curves from S = 0, and add a dashed horizontal line at each Vmax.
3. Round estimates and SE to 2 significant figures for the SE.
4. If both Vmax and Km change, say "consistent with mixed inhibition". If neither changes, say "no clear inhibition detected".
## Correction prompt2
The rounding is wrong. Fix it like this:
- Round the SE to 2 significant figures.
- Round the estimate to the same number of decimal places as the rounded SE.
- Never use scientific notation (no "1.0e2").
Example: Vmax 100.2 with SE 0.9599 should display as "100.20 ± 0.96".
Example: Km 2.008 with SE 0.06037 should display as "2.008 ± 0.060".