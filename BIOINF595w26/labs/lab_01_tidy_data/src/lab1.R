library(tidyverse)

# --- 1. Load and Process Data ---
dr_plate <- readr::read_tsv(
  file = "../data/dr_plate_20210407.txt",
  skip = 34,
  n_max = 58,
  show_col_types = FALSE
) |>
  dplyr::rename(
    time_point = 1,
    temperature = 2
  ) |>
  tidyr::pivot_longer(
    cols = -c(time_point, temperature),
    names_to = "well_id",
    values_to = "growth"
  ) |>
  dplyr::mutate(
    row_letter = stringr::str_extract(well_id, "[A-Z]"),
    column = as.integer(stringr::str_extract(well_id, "[0-9]+")),
    row = match(row_letter, LETTERS),
    
    treatment = dplyr::case_when(
      row %in% c(1, 8) ~ "DMSO",
      row %in% c(2, 3, 4) ~ "BAY 11-7082",
      row %in% c(5, 6, 7) ~ "BAY 11-7085"
    ),
    
    dose_uM = 100 / (2^(column - 1)),
    log10_dose = log10(dose_uM * 1e-6),
    
    dose_label = paste0(signif(dose_uM, 2), " uM"),
    dose_label = factor(
      dose_label,
      levels = paste0(signif(sort(unique(dose_uM), decreasing = TRUE), 2), " uM")
    ),
    # Convert time to numeric hours
    time_point = as.numeric(time_point) / 3600
  )

# Create replica index
replica_map <- dr_plate |>
  dplyr::distinct(row, treatment) |>
  dplyr::group_by(treatment) |>
  dplyr::mutate(replica = factor(dplyr::row_number())) |>
  dplyr::ungroup() |>
  dplyr::select(row, treatment, replica)

# Join replica and select FINAL columns (must keep time_point and replica)
dr_plate <- dr_plate |>
  dplyr::left_join(replica_map, by = c("row", "treatment")) |>
  dplyr::select(
    row, column, treatment, log10_dose, dose_label, 
    temperature, growth, replica, time_point
  )

# Save Data
if (!dir.exists("intermediate")) dir.create("../intermediate")
save(dr_plate, file = "../intermediate/dr_plate.Rdata")

# --- 2. Plot 1: Growth Curves ---
if (!dir.exists("../product")) dir.create("../product")

p1 <- dr_plate |>
  ggplot2::ggplot(ggplot2::aes(x = time_point, y = growth * 100, color = replica)) +
  ggplot2::geom_line() +
  ggplot2::facet_grid(treatment ~ dose_label) +
  ggplot2::scale_x_continuous(breaks = c(0, 6, 12, 18, 24)) +
  ggplot2::scale_y_continuous(breaks = c(0, 50, 100)) +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
  ggplot2::labs(
    x = "Time Point (h)",
    y = "OD of Control (%)",
    color = "Replica",
    title = "Drug Inhibition of Growth Kinetics"
  )

ggplot2::ggsave(
  filename = paste0("../product/growth_curves_", format(Sys.Date(), "%Y%m%d"), ".pdf"),
  plot = p1,
  width = 12, height = 8
)

# --- 3. Plot 2: Dose Response at End-Points ---

# Create a map for axis labels
dose_axis_map <- dr_plate |>
  dplyr::distinct(log10_dose, dose_label) |>
  dplyr::arrange(log10_dose)

p2 <- dr_plate |>
  dplyr::filter(time_point %in% c(6, 12, 18, 24)) |>
  ggplot2::ggplot(ggplot2::aes(x = log10_dose, y = growth * 100, color = replica)) +
  ggplot2::geom_line() +
  ggplot2::facet_grid(time_point ~ treatment) +
  ggplot2::scale_x_continuous(
    breaks = dose_axis_map$log10_dose, 
    labels = dose_axis_map$dose_label
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
  ggplot2::labs(
    x = "Concentration",
    y = "OD of Control (%)",
    color = "Replica",
    title = "Dose Response at Different End-Points"
  )

ggplot2::ggsave(
  filename = paste0("../product/dose_response_", format(Sys.Date(), "%Y%m%d"), ".pdf"),
  plot = p2,
  width = 10, height = 8
)