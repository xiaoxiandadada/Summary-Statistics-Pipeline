plot_manhattan <- function(df, out_path, title = "Manhattan Plot", threshold = NULL,
                           chr_column = "chr", pos_column = "pos", value_column = "W",
                           selected_column = "selected", threshold_column = "W_threshold") {
  required_cols <- c(chr_column, pos_column, value_column)
  if (!all(required_cols %in% names(df))) {
    stop("缺少必要列: ", paste(required_cols, collapse = ", "))
  }

  data <- df[!is.na(df[[chr_column]]) & !is.na(df[[pos_column]]) & !is.na(df[[value_column]]), , drop = FALSE]
  if (nrow(data) == 0) {
    warning("输入数据为空，未生成图像")
    return(invisible())
  }

  data$CHR <- as.integer(data[[chr_column]])
  data$BP <- as.numeric(data[[pos_column]])
  data$VALUE <- as.numeric(data[[value_column]])
  data$SNP <- if ("id" %in% names(data)) as.character(data$id) else paste0(data$CHR, ":", data$BP)
  if (!selected_column %in% names(data)) data[[selected_column]] <- FALSE

  if (is.null(threshold) && threshold_column %in% names(data)) {
    thresholds <- unique(stats::na.omit(as.numeric(data[[threshold_column]])))
    if (length(thresholds) > 0) threshold <- max(thresholds)
  }

  chr_levels <- sort(unique(data$CHR))
  colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b")
  data$COLOR <- colors[(match(data$CHR, chr_levels) - 1) %% length(colors) + 1]

  data$CUM_POS <- NA_real_
  tick_pos <- numeric(length(chr_levels))
  cumulative <- 0
  for (i in seq_along(chr_levels)) {
    chr_val <- chr_levels[i]
    idx <- which(data$CHR == chr_val)
    if (length(idx) == 0) next
    bp <- data$BP[idx]
    data$CUM_POS[idx] <- bp + cumulative
    tick_pos[i] <- cumulative + stats::median(bp)
    cumulative <- cumulative + max(bp)
  }

  y_vals <- data$VALUE
  if (!is.null(threshold) && is.finite(threshold)) y_vals <- c(y_vals, threshold)
  y_lim <- range(y_vals, finite = TRUE)
  y_pad <- diff(y_lim)
  if (!is.finite(y_pad) || y_pad == 0) y_pad <- abs(y_lim[1]) * 0.1 + 1
  y_lim <- c(y_lim[1] - 0.05 * y_pad, y_lim[2] + 0.05 * y_pad)

  png(filename = out_path, width = 1600, height = 600)
  on.exit(dev.off(), add = TRUE)
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  par(mar = c(5, 5, 3, 1))

  plot(
    data$CUM_POS, data$VALUE,
    col = data$COLOR, pch = 16, cex = 1.1,
    xaxt = "n", xlab = "Chromosome", ylab = value_column,
    main = title, cex.axis = 1.2, cex.lab = 1.4, cex.main = 1.6,
    las = 1, font.main = 2, font.lab = 2, ylim = y_lim
  )
  axis(1, at = tick_pos, labels = chr_levels, cex.axis = 1.2)
  axis(2, las = 1, cex.axis = 1.2)
  box()

  if (selected_column %in% names(data)) {
    highlight_idx <- which(data[[selected_column]] %in% TRUE)
    if (length(highlight_idx) > 0) {
      points(
        data$CUM_POS[highlight_idx], data$VALUE[highlight_idx],
        col = "#e41a1c", pch = 16, cex = 1.3
      )
    }
  }

  if (!is.null(threshold) && is.finite(threshold)) {
    abline(h = threshold, col = "#e41a1c", lwd = 2, lty = 2)
  }

  legend_entries <- "All data"
  legend_cols <- "#1f77b4"
  legend_pch <- 16
  legend_lty <- NA
  legend_lwd <- NA

  if (selected_column %in% names(data) && any(data[[selected_column]] %in% TRUE)) {
    legend_entries <- c(legend_entries, "Selected")
    legend_cols <- c(legend_cols, "#e41a1c")
    legend_pch <- c(legend_pch, 16)
    legend_lty <- c(legend_lty, NA)
    legend_lwd <- c(legend_lwd, NA)
  }
  if (!is.null(threshold) && is.finite(threshold)) {
    legend_entries <- c(legend_entries, sprintf("Threshold = %.3f", threshold))
    legend_cols <- c(legend_cols, "#e41a1c")
    legend_pch <- c(legend_pch, NA)
    legend_lty <- c(legend_lty, 2)
    legend_lwd <- c(legend_lwd, 2)
  }

  legend(
    "topright",
    legend = legend_entries,
    col = legend_cols,
    pch = legend_pch,
    lty = legend_lty,
    lwd = legend_lwd,
    bty = "n",
    cex = 1.0
  )

  invisible(out_path)
}
