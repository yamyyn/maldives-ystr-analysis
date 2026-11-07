qc_config <- initial_qc_defaults(single = c("dys576", "dys389i", "dys635", "dys389ii", "dys627", "dys460", "dys458", "dys19", 
"ygatah4", "dys448", "dys391", "dys456", "dys390", "dys438", "dys392", "dys518", "dys570", "dys437", "dys449", "dys393", 
"dys439", "dys481", "dys533"), multicopy = list(dys385_state = c("dys385a", "dys385b"), dyf387s1_state = c("dyf387s1a", "dyf387s1b", 
    "dyf387s1c")))
qc_config$population_mapping <- read_qc_csv("data/example_synthetic/synthetic_population_mapping.csv")
qc_config$confirmed_patterns <- read_qc_csv("data/example_synthetic/synthetic_confirmations.csv")
qc_input <- "data/example_synthetic/initial_qc_synthetic.csv"
qc_output <- "outputs/preparation_example"
