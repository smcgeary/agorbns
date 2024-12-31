library(rje)


ColorViridisPalette <- function(values,
  steps, min, max, log=FALSE, nacolor="gray90", mincolor=NULL, maxcolor=NULL,
  start=0.0, r=0.4, hue=0.8, palettemax=1
) {
  color_ramp <- rev(cubeHelix(round(steps*(1/palettemax)), start=start, r=r, hue=hue))[1:steps]
  if (class(nacolor) == "character") {
    color_ramp <- c(color_ramp, nacolor)
  }
  if (class(mincolor) == "character") {
    color_ramp[1] <- mincolor
  }
  if (class(maxcolor) == "character") {
    color_ramp[steps] <- maxcolor
  }
  if (log) {
    values <- log10(values)
    min <- log10(min)
    max <- log10(max)
  }
  # Transform the data such that the max and min values are 0 and 1
  values <- (values - min)/(max - min)
  # Further transform the data such that the max and min values are `1` and
  # `steps`, respectively.
  col_inds <- round(values*(steps - 1) + 1)
  col_inds <- sapply(col_inds, function(col_ind) {
    min(max(1, col_ind), steps)
  })
  col_inds[which(is.na(col_inds))] <- steps + 1
  return(color_ramp[col_inds])
}