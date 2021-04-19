#' Calculate clonotype abundance
#'
#' @param SCE_in SCE object containing V(D)J data (or data.frame'able object)
#' @param clonotype_col meta.data column containing clonotype IDs
#' @param cluster_col meta.data column containing cluster IDs to use for
#' grouping cells when calculating clonotype abundance
#' @param prefix Prefix to add to new meta.data columns
#' @param return_SCE Return an SCE object. If set to FALSE, a tibble
#' summarizing the results is returned.
#' @return Seurat object with clonotype abundance added to meta.data
#' @import SingleCellExperiment
#' @author djvdj authors
#' @export
calc_abundance <- function(SCE_in, clonotype_col = "cdr3_nt", cluster_col = NULL,
    prefix = "", return_SCE = TRUE) {

    # Format meta.data
    if(any(class(SCE_in) == "SingleCellExperiment")){
        meta_df <- .extract_colData(SCE_in, clonotype_col)
    }else{
        meta_df <- as.data.frame(SCE_in)
        meta_df <- .wrestle_tcr_df(meta_df, clonotype_col)
        return_SCE <- FALSE
    }
    meta_df <- dplyr::select(meta_df, .data$.cell_id, all_of(c(cluster_col, clonotype_col)))

    # Calculate clonotype abundance
    meta_df <- .calc_abund(
        df_in     = meta_df,
        cell_col  = ".cell_id",
        clone_col = clonotype_col,
        clust_col = cluster_col
    )

    # Add results to meta.data
    if (!return_SCE) {
        res <- dplyr::select(meta_df, -.data$.cell_id)
        res <- dplyr::distinct(res)

        return(res)
    }

    new_cols <- c("freq", "pct")

    if (!is.null(cluster_col)) {
        new_cols <- c(new_cols, "shared")
    }

    new_cols <- purrr::set_names(
        new_cols,
        paste0(prefix, "clone_", new_cols)
    )

    meta_df <- select(
        meta_df,
        -all_of(c("n_cells", cluster_col)),
        !!!syms(new_cols)
    )

    meta_df <- tibble::column_to_rownames(meta_df, ".cell_id")

    if(return_SCE==TRUE){
        res <- .add_colData(SCE_in, meta_df)#Seurat::AddMetaData(SCE_in, metadata = meta_df)
    }else{
        res <- meta_df
    }

    return(res)
}


#' Calculate clonotype abundance
#'
#' @param df_in Input data.frame
#' @param cell_col Column containing cell IDs
#' @param clone_col Column containing clonotype IDs
#' @param clust_col Column containing cluster IDs
#' @author djvdj authors
#' @return data.frame
.calc_abund <- function(df_in, cell_col, clone_col, clust_col = NULL) {

    # Count number of cells in each group
    if (!is.null(clust_col)) {
        df_in <- dplyr::group_by(df_in, !!sym(clust_col))
    }

    df_in <- dplyr::mutate(
        df_in,
        n_cells = dplyr::n_distinct(!!sym(cell_col))
    )

    # Calculate frequency
    res <- dplyr::group_by(df_in, !!sym(clone_col), .add = TRUE)

    res <- dplyr::mutate(
        res,
        freq = dplyr::n_distinct(!!sym(cell_col)),
        pct  = (.data$freq / .data$n_cells) * 100
    )

    # Identify shared clonotypes
    if (!is.null(clust_col)) {
        res <- dplyr::group_by(res, !!sym(clone_col))

        res <- dplyr::mutate(
            res,
            shared = dplyr::n_distinct(!!sym(clust_col)) > 1
        )
    }

    res <- dplyr::ungroup(res)

    res
}
