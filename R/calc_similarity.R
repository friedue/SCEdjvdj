#' Calculate repertoire similarity between clusters
#'
#' @param SCE_in Seurat object
#' @param clonotype_col meta.data column containing clonotype IDs to use for
#' calculating overlap
#' @param cluster_col meta.data column containing cluster IDs to use for
#' calculating overlap
#' @param method Method to use for calculating similarity between clusters
#' @param prefix Prefix to add to new meta.data columns
#' @param return_SCE Return a Seurat object. If set to FALSE, a matrix is
#' returned
#' @author djvdj authors
#' @return Seurat object with similarity index added to meta.data
#' @export
calc_similarity <- function(SCE_in, clonotype_col = "cdr3_nt",
    cluster_col, method = abdiv::jaccard,
    prefix = NULL, return_SCE = TRUE) {

    if (is.null(prefix)) {
        prefix <- as.character(substitute(method))
        prefix <- dplyr::last(prefix)
        prefix <- paste0(prefix, "_")
    }

    # Format meta.data
    meta_df <- .extract_colData(SCE_in, clonotype_col)
    meta_df <- dplyr::select(meta_df,
        all_of(c(".cell_id", clonotype_col, cluster_col))
    )

    vdj_df <- dplyr::group_by(
        meta_df,
        !!!syms(c(cluster_col, clonotype_col))
    )

    vdj_df <- dplyr::summarize(
        vdj_df,
        n       = n_distinct(.data$.cell_id),
        .groups = "drop"
    )

    vdj_df <- tidyr::pivot_wider(
        vdj_df,
        names_from  = all_of(cluster_col),
        values_from = .data$n
    )

    # Calculate similarity index
    clusts <- colnames(vdj_df)
    clusts <- clusts[clusts != clonotype_col]

    vdj_df <- dplyr::mutate(vdj_df, across(
        all_of(clusts), ~ tidyr::replace_na(.x, 0)
    ))

    combs <- utils::combn(clusts, 2, simplify = FALSE)

    res <- purrr::map_dfr(combs, ~ {
        ins <- paste0("vdj_df$", .x)

        x <- pull(vdj_df, .x[1])
        y <- pull(vdj_df, .x[2])

        tibble::tibble(
            Var1 = .x[1],
            Var2 = .x[2],
            sim  = method(x, y)
        )
    })

    res_i <- dplyr::rename(res, Var1 = .data$Var2, Var2 = .data$Var1)

    # Return matrix
    if (!return_SCE) {
        res <- tidyr::pivot_wider(
            res,
            names_from  = .data$Var1,
            values_from = .data$sim
        )

        res <- tibble::column_to_rownames(res, "Var2")
        res <- as.matrix(res)

        return(res)
    }

    # Add inverse combinations
    clusts <- tibble::tibble(Var1 = clusts, Var2 = clusts, sim = 1)
    res    <- dplyr::bind_rows(res, res_i, clusts)
    res    <- dplyr::mutate(res, Var1 = paste0(prefix, .data$Var1))

    res <- tidyr::pivot_wider(
        res,
        names_from  = .data$Var1,
        values_from = .data$sim
    )

    # Add similarity index to meta.data
    j_cols <- purrr::set_names("Var2", cluster_col)

    meta_df <- dplyr::left_join(meta_df, res, by = j_cols)
    meta_df <- dplyr::select(meta_df, -all_of(cluster_col))
    meta_df <- tibble::column_to_rownames(meta_df, ".cell_id")
    meta_df <- as.data.frame(meta_df)

    res <- .add_colData(SCE_in, meta_df)
    return(res)
}
