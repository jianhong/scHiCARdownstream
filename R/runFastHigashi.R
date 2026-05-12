#' Run FastHigashi
#' @param inputDirectory The file path for Higashi
#' @param use_conda Use conda to create the python enviroment or not
#' @return nothing.
#' @importFrom reticulate conda_binary install_miniconda import py_save_object
#'  conda_create conda_install conda_list use_condaenv py_config
#' @export
#' @examples
#' # example code
#'
runFastHigashi <- function(inputDirectory='per_cell', use_conda = FALSE){
  # Ensure higashi and fasthigashi are installed in your python env
  # py_install(c("higashi", "fasthigashi", "umap-learn"))

  if(use_conda){
    prepare_miniconda()
    prepare_env()
  }

  higashi <- import("higashi.Higashi_wrapper")
  fasthigashi <- import("fasthigashi.FastHigashi_Wrapper")
  umap <- import("umap")
  pickle <- import("pickle")

  # 2. Process Input Data
  config_path <- file.path(inputDirectory, "Fasthigashi.JSON")
  higashi_model <- higashi$Higashi(config_path)
  higashi_model$process_data()

  # 3. Initialize Fast-Higashi
  config_fast <- file.path(inputDirectory, "Fasthigashi.JSON")
  fastHigashi_model <- fasthigashi$FastHigashi(
    config_path = config_fast,
    path2input_cache = file.path(inputDirectory, "Fasthigashi"),
    path2result_dir = file.path(inputDirectory, "Fasthigashi"),
    off_diag = 100L,
    filter = FALSE,
    do_conv = FALSE,
    do_rwr = FALSE,
    do_col = FALSE,
    no_col = FALSE
  )
  fastHigashi_model$prep_dataset(batch_norm = TRUE)

  # 4. Run the Model
  fastHigashi_model$run_model(
    dim1 = 0.6,
    rank = 256L,
    n_iter_parafac = 1L,
    extra = ""
  )

  # Save model
  py_save_object(fastHigashi_model,
                 file.path(inputDirectory, "fasthigashi_model.pkl"))

  # 5. Run UMAP
  embed_data <- fastHigashi_model$fetch_cell_embedding(final_dim = 256L,
                                                       restore_order = FALSE)
  embedding <- embed_data[['embed_l2_norm_correct_coverage_fh']]

  reducer <- umap$UMAP(n_components = 2L, n_neighbors = 25L, random_state = 0L)
  vec <- reducer$fit_transform(embedding)

  # Convert the matrix to a data frame and add column names
  vec_df <- as.data.frame(vec)
  colnames(vec_df) <- c("UMAP_1", "UMAP_2")

  # Save as a TSV table
  write.table(vec_df,
              file = file.path(inputDirectory, "embedding_coordinates.tsv"),
              sep = "\t",
              row.names = TRUE,
              col.names = TRUE,
              quote = FALSE)

  # 6. Train Higashi for Imputation
  higashi_model_imp <- higashi$Higashi(config_path)
  higashi_model_imp$prep_model()
  higashi_model_imp$train_for_imputation_nbr_0()
  higashi_model_imp$impute_no_nbr()
  higashi_model_imp$train_for_imputation_with_nbr()
  higashi_model_imp$impute_with_nbr()

  # Save model
  py_save_object(higashi_model_imp, file.path(inputDirectory, "higashi_model.impute.pkl"))
}


prepare_miniconda <- function(){
  # 1. Check if Miniconda (or any Conda) is already installed
  is_conda_available <- tryCatch({
    # This returns the path to the conda binary if found
    !is.null(reticulate::conda_binary())
  }, error = function(e) {
    FALSE
  })

  # 2. Install Miniconda if it's missing
  if (!is_conda_available) {
    message("Miniconda not found. Starting installation...")
    reticulate::install_miniconda()
  }

  # 3. Final verification
  py_config()
}

prepare_env <- function(){
  # 1. Define environment name
  env_name <- "higashi_env"

  # 2. Automatically create environment if it doesn't exist
  if (!env_name %in% conda_list()$name) {
    message("Creating conda environment: ", env_name)
    conda_create(env_name, packages = c("python=3.9", "numpy=1.22.3", "pandas=1.3.5", "higashi", "fasthigashi", "umap-learn"), channel = c('ruochiz', 'conda-forge')) # Higashi typically needs Python 3.8+
  }

  # 4. Force R to use this specific environment
  use_condaenv(env_name, required = TRUE)

  # Verify setup
  py_config()
}
