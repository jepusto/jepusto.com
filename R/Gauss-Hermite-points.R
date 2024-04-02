library(rmutil)
library(jsonlite)

lapply(1:10, gauss.hermite) |>
  toJSON() |>
  write("R/Gauss-Hermite-points.json")
