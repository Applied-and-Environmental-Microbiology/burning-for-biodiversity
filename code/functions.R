
### Summarise model coefficients ###############################################
summarise_results <- function(
    model, test_type, exponentiate = FALSE, standardize = NULL
    ) {
  summary_data <- parameters(
    model,
    exponentiate = exponentiate,
    standardize = standardize
  ) %>%
    as_tibble() %>%
    mutate(
      Parameter = ifelse(
        Parameter == "(Intercept)", "(Intercept) [unburnt control]",
        ifelse(grepl("treatment", Parameter),
               sub("^(.*\\btreatment)([A-Za-z])(\\d+)(.*)$",
                   "\\1 [\\2\\3]\\4", Parameter),
               Parameter)),
      across(c(Coefficient, SE, CI, CI_low, CI_high, all_of(test_type)),
             ~round(., 2))
    )
  
  return(summary_data)
}


### Save test statistics #######################################################
save_results_to_excel <- function(response_variables, output_file, overwrite = FALSE) {
  if (file.exists(output_file) && !overwrite) {
    stop("File already exists. To overwrite, set 'overwrite = TRUE'")
  }
  
  # Create a new workbook
  wb <- createWorkbook()
  
  # Loop over response variables
  for (response_var in response_variables) {
    # Object to save
    summary_obj <- get(paste0("summary_", response_var))
    estimate_obj <- get(paste0("estimate_", response_var))
    pairwise_obj <- get(paste0("pairwise_", response_var))
    
    # Add sheets to the workbook for each response variable
    addWorksheet(wb, sheetName = paste0("summary_", response_var))
    writeData(wb, sheet = paste0("summary_", response_var), x = summary_obj)
    
    addWorksheet(wb, sheetName = paste0("estimate_", response_var))
    writeData(wb, sheet = paste0("estimate_", response_var), x = estimate_obj)
    
    addWorksheet(wb, sheetName = paste0("pairwise_", response_var))
    writeData(wb, sheet = paste0("pairwise_", response_var), x = pairwise_obj)
  }
  
  # Save the workbook to a file
  saveWorkbook(wb, output_file, overwrite = TRUE)
}

### My ggplot theme ############################################################
MyTheme = function() {
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.25),
    # The coord_cartesian() function has been implemented to allow annotations
    # outside of the panel boarders. However, that it is boldening the panel
    # boarder, so I have adjusted linewidth here
    axis.text = element_text(colour = 'black', size = rel(0.6)),
    axis.title = element_text(size = rel(0.8)),
    axis.ticks = element_blank(),
    legend.position = "none"
  )
}
