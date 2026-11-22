# Evaluation

## Did the agent give me what I wanted?
Mostly yes, after corrections. The app fits the data, plots it correctly,
and identifies the inhibition type in all the tests I ran.

## Tests

### Test 1: Competitive inhibition (`enzyme_kinetics.csv`)
| Condition | Vmax (µmol/min) | Km (mM) |
|---|---|---|
| control | 100.17 ± 0.96 | 2.008 ± 0.060 |
| inhibitor | 100.8 ± 1.2 | 5.08 ± 0.14 |

- Expected: Vmax unchanged, Km increased → competitive
- App said: "consistent with competitive inhibition"
- PASS
- Only one dashed Vmax line is visible because both Vmax values are about 100
  and the lines overlap. That is the expected result, not a bug.

### Test 2: Non-competitive and mixed (`test_types.csv`, shown as "test2")
| Condition | Vmax | Km | Expected | App said |
|---|---|---|---|---|
| control | 100.17 ± 0.96 | 2.008 ± 0.060 | (reference) | n/a |
| noncomp | 60.29 ± 0.59 | 2.017 ± 0.062 | non-competitive | non-competitive ✅ |
| mixed | 60.35 ± 0.77 | 5.05 ± 0.15 | mixed | mixed ✅ |

- PASS
- Colours match between table and plot for all three conditions.
- The noncomp and mixed Vmax lines overlap at 60, so only one (blue) line shows.

### Test 3: Missing data (`test_missing.csv`)
- Expected: the row with the missing value is dropped and the user is told
- App showed: "1 row dropped" warning, table says "7 read · 6 used",
  plot shows 6 points, fit gives 100.0 ± 1.1 / 1.995 ± 0.071
- With only one condition, the app explains that it needs two conditions
  to compare, instead of crashing.
- PASS

## What worked
- Nonlinear fit with `nls()`, standard errors reported
- Correct inhibition type in all three tested cases
- Consistent colours, y-axis from 0, curves from the origin, Vmax lines
- Bad rows handled safely with a warning
- Runs fully in the browser (webR), so no data is uploaded

## What did not work / limitations
- First rounding attempt was wrong (`1.0e2`); fixed with a clearer prompt
- The dropped-row warning doesn't say which row was dropped or why
- Overlapping Vmax lines can hide one condition's line
- Not tested:
uncompetitive inhibition, the "no inhibition" message,
  wrong column names, a fit that fails to converge
- "Similar" means within one combined SE. This is a simple rule of thumb,
  not a statistical test, so noisy real data could be misclassified.
- The first condition is always treated as the control
- No Ki is calculated; only a qualitative inhibition type is given

## Evidence
- R code panel showing `nls()`: ![R code]!![alt text](<螢幕截圖 2026-11-22 下午12.41.53.png>)

## Full R code
input <- read.csv("/input.csv", stringsAsFactors = FALSE,
                  na.strings = c("", "NA"), check.names = FALSE)
required <- c("condition", "substrate_mM", "velocity")
if (!all(required %in% names(input))) {
  stop(paste("CSV must contain:", paste(required, collapse = ", ")))
}

input$condition <- trimws(as.character(input$condition))
input$substrate_mM <- suppressWarnings(as.numeric(input$substrate_mM))
input$velocity <- suppressWarnings(as.numeric(input$velocity))
condition_label <- input$condition
condition_label[is.na(condition_label) | condition_label == ""] <- "(missing condition)"
condition_order <- unique(condition_label)
condition_colors <- grDevices::hcl.colors(length(condition_order), palette = "Dark 3")
row_counts <- data.frame(
  condition = condition_order,
  rows_read = vapply(condition_order, function(x) sum(condition_label == x), integer(1)),
  color = condition_colors
)

valid <- !is.na(input$condition) & input$condition != "" &
  is.finite(input$substrate_mM) & input$substrate_mM >= 0 &
  is.finite(input$velocity)
data <- input[valid, required]
write.table(row_counts, "/counts.tsv", sep = intToUtf8(9), row.names = FALSE, quote = TRUE)

fits <- lapply(unique(data$condition), function(condition) {
  d <- data[data$condition == condition, , drop = FALSE]
  tryCatch({
    fit <- nls(
      velocity ~ Vmax * substrate_mM / (Km + substrate_mM),
      data = d,
      start = list(Vmax = max(d$velocity), Km = median(d$substrate_mM))
    )
    estimate <- summary(fit)$parameters
    data.frame(
      condition = condition,
      Vmax = unname(coef(fit)[["Vmax"]]),
      Vmax_SE = estimate["Vmax", "Std. Error"],
      Km = unname(coef(fit)[["Km"]]),
      Km_SE = estimate["Km", "Std. Error"],
      rows_used = nrow(d),
      status = "ok",
      error = ""
    )
  }, error = function(e) {
    data.frame(
      condition = condition, Vmax = NA_real_, Vmax_SE = NA_real_,
      Km = NA_real_, Km_SE = NA_real_, rows_used = nrow(d),
      status = "failed", error = conditionMessage(e)
    )
  })
})
if (length(fits)) {
  results <- do.call(rbind, fits)
} else {
  results <- data.frame(
    condition = character(), Vmax = numeric(), Vmax_SE = numeric(),
    Km = numeric(), Km_SE = numeric(), rows_used = integer(),
    status = character(), error = character()
  )
}
write.table(results, "/results.tsv", sep = intToUtf8(9), row.names = FALSE, quote = TRUE, na = "")

svg("/plot.svg", width = 8, height = 5, pointsize = 11)
if (nrow(data) > 0) {
  conditions <- unique(data$condition)
  colors <- condition_colors[match(conditions, condition_order)]
  fitted_vmax <- results$Vmax[results$status == "ok"]
  y_top <- max(c(0, data$velocity, fitted_vmax), na.rm = TRUE) * 1.1
  if (!is.finite(y_top) || y_top <= 0) y_top <- 1
  x_top <- max(data$substrate_mM)
  if (x_top <= 0) x_top <- 1
  plot(
    data$substrate_mM, data$velocity,
    col = colors[match(data$condition, conditions)], pch = 19,
    xlim = c(0, x_top),
    ylim = c(0, y_top),
    xlab = "Substrate concentration (mM)",
    ylab = "Velocity (µmol/min)",
    main = "Michaelis–Menten kinetics"
  )
  for (i in seq_along(conditions)) {
    fit_row <- match(conditions[i], results$condition)
    if (!is.na(fit_row) && results$status[fit_row] == "ok") {
      curve(
        results$Vmax[fit_row] * x / (results$Km[fit_row] + x),
        from = 0,
        to = max(data$substrate_mM[data$condition == conditions[i]]),
        add = TRUE, col = colors[i], lwd = 2
      )
      abline(h = results$Vmax[fit_row], col = colors[i], lty = 2, lwd = 1)
    }
  }
  fitted <- conditions %in% results$condition[results$status == "ok"]
  legend(
    "bottomright", legend = conditions, col = colors, pch = 19,
    lty = ifelse(fitted, 1, NA), lwd = 2, bty = "n"
  )
} else {
  plot.new()
  title("No valid measurements to plot")
}
dev.off()
