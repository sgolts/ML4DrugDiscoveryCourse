library(httr2)
library(purrr)
library(dplyr)
library(readr)

source("pubchem_utils.R") 

# AIDs - no 592 because no active compounds!
aids <- c(587, 588, 590, 591, 593, 594)

# split apply combine
all_active_compounds <- aids |> 
  purrr::map_dfr(function(current_aid) {
    
    message(paste("Active compounds for AID:", current_aid))
    
    compounds <- get_aid_compounds(
      aid = current_aid,
      cids_type = "active", 
      properties = c("SMILES", "InChI", "InChIKey"),
      verbose = FALSE
    )
    
    compounds |> 
      dplyr::mutate(AID = current_aid)
    
  })

# count actives in multiple datasets
multi_active_counts <- all_active_compounds |>
  dplyr::group_by(CID) |>
  dplyr::summarise(
    n_assays = dplyr::n_distinct(AID),
    assays = paste(unique(AID), collapse = ", ")
  ) |>
  dplyr::filter(n_assays > 1) |>
  dplyr::arrange(desc(n_assays))

# print outputs
cat("\nAnalysis Results\n")
cat("Total active entries retrieved:", nrow(all_active_compounds), "\n")
cat("Number of unique active compounds:", length(unique(all_active_compounds$CID)), "\n")
cat("Number of compounds active in more than 1 dataset:", nrow(multi_active_counts), "\n")

if(nrow(multi_active_counts) > 0) {
  cat("\nCompounds active in most assays:\n")
  print(head(multi_active_counts))
}

# save data
# order: AID, CID, SMILES, InChI, InChIKey
final_output <- all_active_compounds |>
  dplyr::select(AID, CID, SMILES, InChI, InChIKey)

readr::write_tsv(final_output, "../intermediate/active_compounds.tsv")

message("\nData saved to 'intermediate/active_compounds.tsv'")