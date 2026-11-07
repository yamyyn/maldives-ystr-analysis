# Original mathematical/parsing helpers; top-level study scripts were not executed.
run_xval <- function(matrix_data, grouping_factor, n_da, analysis_name, seed) {
    log_both(paste("Starting", analysis_name, "xvalDapc:", length(pc_candidates), "PC candidates x", N_REPETITIONS, "repetitions."))
    set.seed(seed)
    result <- adegenet::xvalDapc(matrix_data, grouping_factor, n.pca.max = matrix_rank, n.da = n_da, training.set = TRAINING_PROPORTION, 
        result = "groupMean", center = TRUE, scale = FALSE, n.pca = pc_candidates, n.rep = N_REPETITIONS, xval.plot = FALSE)
    log_both(paste(analysis_name, "xvalDapc completed."))
    result
}
parse_xval <- function(xval_object, analysis_name) {
    required_names <- c("Cross-Validation Results", "Median and Confidence Interval for Random Chance", "Mean Successful Assignment by Number of PCs of PCA", 
        "Number of PCs Achieving Highest Mean Success", "Root Mean Squared Error by Number of PCs of PCA", "Number of PCs Achieving Lowest MSE", 
        "DAPC")
    if (!all(required_names %in% names(xval_object))) {
        stop(paste(analysis_name, "xvalDapc output is missing expected named components."), call. = FALSE)
    }
    cv_raw <- xval_object[["Cross-Validation Results"]] %>% as_tibble() %>% mutate(Analysis = analysis_name, n.pca = as.integer(n.pca), 
        success = as.numeric(success)) %>% select(Analysis, n.pca, success)
    mean_success <- xval_object[["Mean Successful Assignment by Number of PCs of PCA"]]
    rmse <- xval_object[["Root Mean Squared Error by Number of PCs of PCA"]]
    summary_table <- full_join(tibble(n.pca = as.integer(names(mean_success)), Mean_Group_Assignment_Success = as.numeric(mean_success)), 
        tibble(n.pca = as.integer(names(rmse)), RMSE = as.numeric(rmse)), by = "n.pca") %>% arrange(n.pca) %>% mutate(Analysis = analysis_name, 
        Mean_Group_Assignment_Percent = 100 * Mean_Group_Assignment_Success) %>% select(Analysis, n.pca, Mean_Group_Assignment_Success, 
        Mean_Group_Assignment_Percent, RMSE)
    best_mse_pc <- as.integer(xval_object[["Number of PCs Achieving Lowest MSE"]])
    highest_success_pc <- as.integer(xval_object[["Number of PCs Achieving Highest Mean Success"]])
    random_chance <- xval_object[["Median and Confidence Interval for Random Chance"]]
    list(Raw = cv_raw, Summary = summary_table, Best_MSE_PC = best_mse_pc, Highest_Success_PC = highest_success_pc, Random_Chance = random_chance, 
        DAPC_Object_From_Xval = xval_object[["DAPC"]])
}
validate_cv_counts <- function(parsed, label) {
    repetition_counts <- parsed$Raw %>% count(n.pca, name = "Repetitions")
    if (any(repetition_counts$Repetitions != N_REPETITIONS)) {
        fail_with_log(paste(label, "does not contain exactly", N_REPETITIONS, "replicates for every PC candidate."))
    }
    if (!parsed$Best_MSE_PC %in% parsed$Summary$n.pca) {
        fail_with_log(paste(label, "best-MSE PC value is not present in the tested grid."))
    }
    repetition_counts
}
get_row_for_pc <- function(summary_table, pc) {
    summary_table %>% filter(n.pca == pc)
}
