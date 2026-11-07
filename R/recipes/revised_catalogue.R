# Adapted from the user's original revised_catalogue final builder.
# Adaptations: configurable labels, font, palettes, axes, cohort-neutral sheets and matrix boundaries; see documentation.
write_table <- function(df, filename) utils::write.csv(df, file.path(tab_dir, filename), row.names = FALSE, na = "")
