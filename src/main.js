import { WebR } from "webr";
import "./style.css";

const R_CODE = `input <- read.csv("/input.csv", stringsAsFactors = FALSE,
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
`;

const MAX_FILE_SIZE = 5 * 1024 * 1024;
const webR = new WebR();
let webRReady;
let csvText = "";
let currentFileName = "Loading example data…";
let plotUrl;

const fileInput = document.querySelector("#csv-file");
const analyzeButton = document.querySelector("#analyze-button");
const fileName = document.querySelector("#file-name");
const fileStatus = document.querySelector("#file-status");
const emptyState = document.querySelector("#empty-state");
const results = document.querySelector("#results");
const tableBody = document.querySelector("#results-table");
const warningBanner = document.querySelector("#warning-banner");
const errorBanner = document.querySelector("#error-banner");
const plotImage = document.querySelector("#plot-image");
const plotSection = document.querySelector("#plot-section");
const interpretationText = document.querySelector("#interpretation-text");
const fitBadge = document.querySelector("#fit-badge");

function initWebR() {
  webRReady ??= webR.init();
  return webRReady;
}

function setFile(text, name) {
  csvText = text;
  currentFileName = name;
  fileName.textContent = name;
  fileStatus.classList.remove("file-error");
  analyzeButton.disabled = false;
  errorBanner.hidden = true;
}

function parseDelimited(text, delimiter = "\t") {
  const rows = [];
  let row = [];
  let field = "";
  let quoted = false;
  for (let i = 0; i < text.length; i += 1) {
    const char = text[i];
    if (quoted) {
      if (char === '"' && text[i + 1] === '"') {
        field += '"';
        i += 1;
      } else if (char === '"') {
        quoted = false;
      } else {
        field += char;
      }
    } else if (char === '"' && field === "") {
      quoted = true;
    } else if (char === delimiter) {
      row.push(field);
      field = "";
    } else if (char === "\n" || char === "\r") {
      if (char === "\r" && text[i + 1] === "\n") i += 1;
      row.push(field);
      if (row.some((cell) => cell !== "")) rows.push(row);
      row = [];
      field = "";
    } else {
      field += char;
    }
  }
  if (field !== "" || row.length > 0) {
    row.push(field);
    if (row.some((cell) => cell !== "")) rows.push(row);
  }
  if (rows.length === 0) return [];
  const headers = rows[0];
  return rows.slice(1).map((cells) =>
    Object.fromEntries(headers.map((header, index) => [header, cells[index] ?? ""])),
  );
}

function formatParameter(estimate, standardError) {
  const estimateValue = Number(estimate);
  const standardErrorValue = Number(standardError);
  if (!Number.isFinite(estimateValue) || !Number.isFinite(standardErrorValue)) return "—";

  const roundedStandardError = Number(standardErrorValue.toPrecision(2));
  const decimalPlaces = roundedStandardError === 0
    ? 0
    : Math.max(0, 1 - Math.floor(Math.log10(Math.abs(roundedStandardError))));
  const formatter = new Intl.NumberFormat("en-US", {
    useGrouping: false,
    minimumFractionDigits: decimalPlaces,
    maximumFractionDigits: decimalPlaces,
  });
  return `${formatter.format(estimateValue)} <span class="pm">±</span> ${formatter.format(roundedStandardError)}`;
}

function escapeHtml(value) {
  return String(value).replace(/[&<>"']/g, (char) => ({
    "&": "&amp;",
    "<": "&lt;",
    ">": "&gt;",
    '"': "&quot;",
    "'": "&#39;",
  })[char]);
}

function combinedSe(a, b) {
  return Math.hypot(Number(a), Number(b));
}

function getInterpretation(fitRows) {
  const successful = fitRows.filter((row) => row.status === "ok");
  if (successful.length < 2) {
    return "Fit at least two conditions successfully to compare inhibitor behavior.";
  }
  const reference = successful.find((row) => row.condition.toLowerCase() === "control") ?? successful[0];
  const comparisons = successful.filter((row) => row !== reference);
  return comparisons.map((condition) => {
    const vmSimilar = Math.abs(Number(reference.Vmax) - Number(condition.Vmax)) <=
      combinedSe(reference.Vmax_SE, condition.Vmax_SE);
    const kmSimilar = Math.abs(Number(reference.Km) - Number(condition.Km)) <=
      combinedSe(reference.Km_SE, condition.Km_SE);
    const vmChanged = !vmSimilar;
    const kmChanged = !kmSimilar;
    const kmIncreased = Number(condition.Km) - Number(reference.Km) >
      combinedSe(reference.Km_SE, condition.Km_SE);
    const vmDecreased = Number(reference.Vmax) - Number(condition.Vmax) >
      combinedSe(reference.Vmax_SE, condition.Vmax_SE);
    let result = "does not clearly match either inhibition pattern";
    if (vmChanged && kmChanged) result = "consistent with mixed inhibition";
    else if (vmSimilar && kmIncreased) result = "consistent with competitive inhibition";
    else if (vmDecreased && kmSimilar) result = "consistent with non-competitive inhibition";
    else if (!vmChanged && !kmChanged) result = "no clear inhibition detected";
    return `${condition.condition} vs ${reference.condition}: ${result}.`;
  }).join(" ");
}

function showResults(countRows, fitRows, droppedRows, plotSvg) {
  emptyState.hidden = true;
  results.hidden = false;
  errorBanner.hidden = true;

  const counts = new Map(countRows.map((row) => [
    row.condition,
    { rows: Number(row.rows_read), color: row.color },
  ]));
  const fits = new Map(fitRows.map((row) => [row.condition, row]));
  const allConditions = [...new Set([...counts.keys(), ...fits.keys()])];
  tableBody.innerHTML = allConditions.map((condition) => {
    const fit = fits.get(condition);
    const count = counts.get(condition);
    const color = /^#[\da-f]{6}$/i.test(count?.color ?? "") ? count.color : "#58a182";
    const status = fit?.status === "failed"
      ? `<span class="fit-failed" title="${escapeHtml(fit.error)}">Fit failed</span>`
      : "";
    return `<tr>
      <th scope="row"><span class="condition-dot" style="--condition-color: ${color}"></span>${escapeHtml(condition)} ${status}</th>
      <td>${fit?.status === "ok" ? formatParameter(fit.Vmax, fit.Vmax_SE) : "—"}</td>
      <td>${fit?.status === "ok" ? formatParameter(fit.Km, fit.Km_SE) : "—"}</td>
      <td>${count?.rows ?? 0} read <span class="rows-used">· ${fit?.rows_used ?? 0} used</span></td>
    </tr>`;
  }).join("");

  const failedFits = fitRows.filter((row) => row.status === "failed");
  warningBanner.hidden = droppedRows === 0 && failedFits.length === 0;
  warningBanner.innerHTML = [
    droppedRows > 0 ? `<strong>${droppedRows} row${droppedRows === 1 ? "" : "s"} dropped.</strong> Rows with a missing condition, non-numeric or non-finite value, or negative substrate concentration were excluded.` : "",
    ...failedFits.map((row) => `<strong>${escapeHtml(row.condition)} fit failed:</strong> ${escapeHtml(row.error)}`),
  ].filter(Boolean).join("<br>");

  fitBadge.textContent = `${fitRows.filter((row) => row.status === "ok").length} / ${allConditions.length} fits`;
  interpretationText.textContent = getInterpretation(fitRows);
  plotSection.hidden = !plotSvg;
  if (plotSvg) {
    if (plotUrl) URL.revokeObjectURL(plotUrl);
    plotUrl = URL.createObjectURL(new Blob([plotSvg], { type: "image/svg+xml" }));
    plotImage.src = plotUrl;
  }
}

async function analyze() {
  if (!csvText) return;
  analyzeButton.disabled = true;
  analyzeButton.classList.add("loading");
  analyzeButton.querySelector("span").textContent = "Starting R…";
  warningBanner.hidden = true;
  errorBanner.hidden = true;
  results.hidden = false;
  emptyState.hidden = true;
  tableBody.innerHTML = "";
  plotSection.hidden = true;

  try {
    await initWebR();
    analyzeButton.querySelector("span").textContent = "Fitting curves…";
    await webR.FS.writeFile("/input.csv", new TextEncoder().encode(csvText));
    await webR.evalRVoid(R_CODE);

    const [countsBytes, fitsBytes, plotBytes] = await Promise.all([
      webR.FS.readFile("/counts.tsv"),
      webR.FS.readFile("/results.tsv"),
      webR.FS.readFile("/plot.svg"),
    ]);
    const decoder = new TextDecoder();
    const countsText = decoder.decode(countsBytes);
    const fitsText = decoder.decode(fitsBytes);
    const plotText = decoder.decode(plotBytes);
    const countRows = parseDelimited(countsText);
    const fitRows = parseDelimited(fitsText);
    const droppedRows = countRows.reduce((total, row) => {
      const fit = fitRows.find((candidate) => candidate.condition === row.condition);
      return total + Math.max(0, Number(row.rows_read) - Number(fit?.rows_used ?? 0));
    }, 0);
    showResults(countRows, fitRows, droppedRows, plotText);
  } catch (error) {
    errorBanner.textContent = error instanceof Error ? error.message : String(error);
    errorBanner.hidden = false;
    warningBanner.hidden = true;
    plotSection.hidden = true;
    interpretationText.textContent = "Analysis could not be completed. Check the CSV format and R error above.";
  } finally {
    analyzeButton.disabled = false;
    analyzeButton.classList.remove("loading");
    analyzeButton.querySelector("span").textContent = "Run analysis";
  }
}

fileInput.addEventListener("change", async (event) => {
  const [file] = event.target.files;
  if (!file) return;
  if (file.size > MAX_FILE_SIZE) {
    fileStatus.classList.add("file-error");
    fileName.textContent = "File exceeds the 5 MB limit.";
    analyzeButton.disabled = true;
    return;
  }
  setFile(await file.text(), file.name);
});

analyzeButton.addEventListener("click", analyze);

const uploadZone = document.querySelector(".upload-zone");
uploadZone.addEventListener("dragover", (event) => {
  event.preventDefault();
  uploadZone.classList.add("dragging");
});
uploadZone.addEventListener("dragleave", () => uploadZone.classList.remove("dragging"));
uploadZone.addEventListener("drop", async (event) => {
  event.preventDefault();
  uploadZone.classList.remove("dragging");
  const [file] = event.dataTransfer.files;
  if (!file) return;
  if (file.size > MAX_FILE_SIZE) {
    fileStatus.classList.add("file-error");
    fileName.textContent = "File exceeds the 5 MB limit.";
    analyzeButton.disabled = true;
    return;
  }
  setFile(await file.text(), file.name);
});

async function loadExample() {
  try {
    const response = await fetch("/data/enzyme_kinetics.csv");
    if (!response.ok) throw new Error(`Could not load example data (${response.status}).`);
    setFile(await response.text(), "enzyme_kinetics.csv");
    await analyze();
  } catch (error) {
    fileStatus.classList.add("file-error");
    fileName.textContent = "Example data could not be loaded.";
    errorBanner.textContent = error instanceof Error ? error.message : String(error);
    errorBanner.hidden = false;
    emptyState.hidden = true;
    results.hidden = false;
  }
}

loadExample();
