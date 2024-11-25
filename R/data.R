#' @title Sample Expression Profile
#'
#' @description The sample-wise cell-type proportion matrix is simulated by the Dirichlet distribution function. We set α (the parameter controlling the shape of the distribution) at 3 for all cell types in our simulation to make sure the proportions among samples are basically concentrated towards the center and mutually symmetric. Here we provide an example A matrix for users’ reference, with 30 samples and 5 cell types.
#'
#' @format A data frame with 30 rows and 5 columns
"sample_expression_profile"

#' @title Sample Proportion Matrix
#'
#' @description The cell type specific expression profiles are extracted from GSE73721, an expression profiling on RNA-Seq of human brain cells. Five cell types were extracted from the original data: HepaCAM positive cells (Astrocyte), GalC positive cells (Oligodendrocyte), CD45 positive cells (Myeloid), BSL bound cells (Endothelial), and O4 positive cells (Oligodendrocyte), all from adult temporal lobe. Here we provided an pre-processed example cell-type specific expression profile for users to use as a reference.
#'
#' @format A data frame with 5 rows and 15,041 columns
"sample_proportion_matrix"