#' Filter V(D)J meta.data
#'
#' @param sobj_in Seurat object containing V(D)J data
#' @param filt Condition to use for filtering meta.data
#' @param clonotype_col meta.data column containing clonotype IDs. This column
#' is used to determine which cells have V(D)J data. If clonotype_col is set to
#' NULL, filtering is performed regardless of whether V(D)J data is present for
#' the cell
#' @param filter_cells Should cells be removed from object? If set FALSE
#' (default) V(D)J data will be removed from meta.data but no cells will be
#' removed from the object
#' @param sep Separator to use for expanding meta.data columns
#' @param vdj_cols meta.data columns containing VDJ data to use for filtering.
#' If set to NULL (recommended) columns are automatically selected.
#' @return Seurat object
#' @export
filter_vdj <- function(sobj_in, filt, clonotype_col = "cdr3_nt", filter_cells = FALSE,
                       sep = ";", vdj_cols = NULL) {

  if (!filter_cells && is.null(clonotype_col)) {
    stop("clonotype_col must be provided when filter_cells is set to FALSE.")
  }

  # Identify columns with VDJ data
  meta_df <- sobj_in@meta.data
  meta_df <- tibble::as_tibble(meta_df, rownames = ".cell_id")

  col_list <- .get_vdj_cols(
    df_in     = meta_df,
    clone_col = clonotype_col,
    cols_in   = vdj_cols,
    sep       = sep
  )

  vdj_cols <- col_list$vdj
  sep_cols <- col_list$sep

  # Create list-cols for VDJ columns that contain sep
  if (!purrr::is_empty(sep_cols)) {
    sep_cols <- set_names(
      x  = sep_cols,
      nm = paste0(".", sep_cols)
    )

    meta_df <- .split_vdj(
      df_in    = meta_df,
      sep_cols = sep_cols,
      sep      = sep
    )

    meta_df <- dplyr::rowwise(meta_df)
  }

  # Store results for filtering
  meta_df <- dplyr::mutate(meta_df, .KEEP = {{filt}})
  meta_df <- dplyr::ungroup(meta_df)

  # Remove VDJ data from meta.data (without filtering cells)
  if (!filter_cells) {
    meta_df <- dplyr::mutate(
      meta_df,
      across(
        all_of(c(vdj_cols, names(sep_cols))),
        ~ ifelse(.KEEP, .x, NA)
      )
    )

  # Filter cells from object
  # In clonotype_col != NULL only filter cells with VDJ data
  } else {
    if (!is.null(clonotype_col)) {
      meta_df <- dplyr::mutate(
        meta_df,
        .KEEP = dplyr::if_else(
          is.na(!!sym(clonotype_col)),
          TRUE,
          .data$.KEEP
        )
      )
    }

    meta_df <- dplyr::filter(meta_df, .data$.KEEP)
  }

  # Remove columns created for filtering
  if (!purrr::is_empty(sep_cols)) {
    sep_cols <- purrr::set_names(
      x  = names(sep_cols),
      nm = unname(sep_cols)
    )

    meta_df <- dplyr::select(
      meta_df,
      !all_of(c(names(sep_cols), ".KEEP"))
    )

    meta_df <- dplyr::rename(meta_df, !!!syms(sep_cols))
  }

  # Add meta.data to Seurat object
  meta_df <- tibble::column_to_rownames(meta_df, ".cell_id")

  cells <- rownames(meta_df)
  res   <- subset(sobj_in, cells = cells)
  res   <- Seurat::AddMetaData(res, meta_df)

  res
}


#' Mutate V(D)J meta.data
#'
#' @param sobj_in Seurat object containing V(D)J data
#' @param ... Name-value pairs to use for creating or modifying meta.data
#' columns
#' @param clonotype_col meta.data column containing clonotype IDs. This column
#' is used to determine which cells have V(D)J data
#' @param sep Separator to use for expanding meta.data columns
#' @param vdj_cols meta.data columns containing VDJ data to use for filtering.
#' If set to NULL (recommended) columns are automatically selected.
#' @return Seurat object
#' @export
mutate_vdj <- function(sobj_in, ..., clonotype_col = "cdr3_nt", sep = ";", vdj_cols = NULL) {

  # Identify columns with VDJ data
  meta_df <- sobj_in@meta.data
  meta_df <- tibble::as_tibble(meta_df, rownames = ".cell_id")

  col_list <- .get_vdj_cols(
    df_in     = meta_df,
    clone_col = clonotype_col,
    cols_in   = vdj_cols,
    sep       = sep
  )

  vdj_cols <- col_list$vdj
  sep_cols <- col_list$sep

  # Create list-cols for VDJ columns that contain sep
  if (!purrr::is_empty(sep_cols)) {
    sep_cols <- set_names(
      x  = sep_cols,
      nm = paste0(".", sep_cols)
    )

    meta_df <- .split_vdj(
      df_in    = meta_df,
      sep_cols = sep_cols,
      sep      = sep
    )

    meta_df <- dplyr::rowwise(meta_df)
  }

  # Mutate meta.data
  meta_df <- dplyr::mutate(meta_df, ...)
  meta_df <- dplyr::ungroup(meta_df)

  # Remove columns created for mutate
  if (!purrr::is_empty(sep_cols)) {
    sep_cols <- purrr::set_names(
      x  = names(sep_cols),
      nm = unname(sep_cols)
    )

    meta_df <- dplyr::select(
      meta_df,
      !all_of(names(sep_cols))
    )

    meta_df <- dplyr::rename(meta_df, !!!syms(sep_cols))
  }

  # Add meta.data to Seurat object
  meta_df <- tibble::column_to_rownames(meta_df, ".cell_id")

  cells <- rownames(meta_df)
  res   <- subset(sobj_in, cells = cells)
  res   <- Seurat::AddMetaData(res, meta_df)

  res
}



#' Calculate repertoire diversity
#'
#' @param sobj_in Seurat object
#' @param clonotype_col meta.data column containing clonotype IDs to use for
#' calculating diversity
#' @param cluster_col meta.data column containing cluster IDs to use for
#' calculating diversity. If cluster_col is omitted, diversity index will be
#' calculated for all clonotypes.
#' @param method Method to use for calculating diversity. Can pass also pass
#' a named list to use multiple methods.
#' @param prefix Prefix to add to new meta.data columns
#' @param return_seurat Return a Seurat object. If set to FALSE, the meta.data
#' table is returned.
#' @return Seurat object with diversity index added to meta.data
#' @export
calc_diversity <- function(sobj_in, clonotype_col = "cdr3_nt", cluster_col = NULL,
                           method = abdiv::simpson, prefix = "", return_seurat = TRUE) {

  if (length(method) > 1 && is.null(names(method))) {
    stop("Must include names if using a list of methods.")
  }

  if (length(method) == 1 && is.null(names(method))) {
    nm <- as.character(substitute(method))
    nm <- last(nm)

    method <- purrr::set_names(list(method), nm)
  }

  # Format meta.data
  vdj_cols <- clonotype_col
  meta_df  <- tibble::as_tibble(sobj_in@meta.data, rownames = ".cell_id")
  vdj_df   <- dplyr::filter(meta_df, !is.na(!!sym(clonotype_col)))

  # Count clonotypes
  if (!is.null(cluster_col)) {
    vdj_cols <- c(cluster_col, vdj_cols)
    vdj_df   <- dplyr::group_by(vdj_df, !!sym(cluster_col))
  }

  vdj_df <- dplyr::group_by(vdj_df, !!sym(clonotype_col), .add = TRUE)

  vdj_df <- dplyr::summarize(
    vdj_df,
    .n      = dplyr::n_distinct(.data$.cell_id),
    .groups = "drop"
  )

  # Calculate diversity
  if (!is.null(cluster_col)) {
    vdj_df <- dplyr::group_by(vdj_df, !!sym(cluster_col))
  }

  if (length(method) == 1 && is.null(names(method))) {
    nm <- as.character(substitute(method))
    nm <- dplyr::last(nm)

    method        <- list(method)
    names(method) <- nm
  }

  div_cols      <- paste0(prefix, names(method))
  names(method) <- div_cols

  vdj_df <- purrr::imap_dfr(method, ~ {
    dplyr::mutate(
      vdj_df,
      diversity = .x(.data$.n),
      met       = .y
    )
  })

  vdj_df <- tidyr::pivot_wider(
    vdj_df,
    names_from  = "met",
    values_from = "diversity"
  )

  vdj_df <- dplyr::ungroup(vdj_df)
  vdj_df <- dplyr::select(vdj_df, all_of(c(vdj_cols, div_cols)))

  # Return data.frame
  if (!return_seurat) {
    res <- dplyr::select(vdj_df, all_of(c(cluster_col, div_cols)))
    res <- dplyr::distinct(res)

    if (is.null(cluster_col)) {
      res <- unlist(res)
    }

    return(res)
  }

  # Add results to meta.data
  meta_df <- dplyr::left_join(meta_df, vdj_df, by = vdj_cols)
  meta_df <- tibble::column_to_rownames(meta_df, ".cell_id")
  meta_df <- as.data.frame(meta_df)

  res <- Seurat::AddMetaData(sobj_in, metadata = meta_df)

  res
}



#' Calculate gene usage
#'
#' Gene usage is calculated as the percent of cells that express each V(D)J
#' gene present in gene_cols. Cells that lack V(D)J data and have NAs present
#' in gene_cols are excluded from this calculation.
#'
#' @param sobj_in Seurat object containing V(D)J data
#' @param gene_cols meta.data column containing genes used for each clonotype
#' @param cluster_col meta.data column containing cell clusters to use when
#' calculating gene usage
#' @param chain Chain to use for calculating gene usage. Set to NULL to include
#' all chains.
#' @param chain_col meta.data column containing chains for each cell
#' @param sep Separator to use for expanding gene_cols
#' @return data.frame containing gene usage summary
#' @export
calc_usage <- function(sobj_in, gene_cols, cluster_col = NULL, chain = NULL,
                       chain_col = "chains", sep = ";") {

  # data.frame to calculate usage
  sep_cols <- c(gene_cols, chain_col)
  vdj_cols   <- c(".cell_id", cluster_col, sep_cols)

  meta_df <- sobj_in@meta.data
  meta_df <- tibble::as_tibble(meta_df, rownames = ".cell_id")
  meta_df <- dplyr::select(meta_df, all_of(vdj_cols))

  meta_df <- dplyr::filter(meta_df, across(
    all_of(gene_cols),
    ~ !is.na(.x)
  ))

  res <- dplyr::mutate(meta_df, across(
    all_of(sep_cols),
    ~ strsplit(as.character(.x), sep)
  ))

  # Filter chains
  if (!is.null(chain)) {
    if (is.null(chain_col)) {
      stop("Must specify chain_col.")
    }

    res <- dplyr::mutate(res, across(all_of(gene_cols), ~ {
      g_col <- .x

      purrr::map2(g_col, !!sym(chain_col), ~ {
        .x <- dplyr::if_else(
          any(.y %in% chain),
          list(.x[.y %in% chain]),
          list("None")
        )

        unlist(.x)
      })
    }))

    res <- dplyr::select(res, -all_of(chain_col))
  }

  res <- tidyr::unnest(res, cols = all_of(gene_cols))
  res <- dplyr::distinct(res)

  # Count genes used
  res <- dplyr::group_by(res, !!!syms(gene_cols))

  if (!is.null(cluster_col)) {
    res <- dplyr::group_by(res, !!sym(cluster_col), .add = TRUE)
  }

  res <- dplyr::summarize(
    res,
    n_cells = n_distinct(meta_df$.cell_id),
    freq    = dplyr::n_distinct(.data$.cell_id),
    .groups = "drop"
  )

  # Calculate percentage used
  if (!is.null(cluster_col)) {
    clusts       <- pull(meta_df, cluster_col)
    clust_counts <- table(clusts)
    clusts       <- unique(clusts)

    res <- tidyr::pivot_wider(
      res,
      names_from  = all_of(cluster_col),
      values_from = .data$freq,
      values_fill = 0
    )

    res <- tidyr::pivot_longer(
      res,
      cols      = all_of(clusts),
      names_to  = cluster_col,
      values_to = "freq"
    )

    res <- dplyr::mutate(  # why doesn't .before work here??
      res,
      n_cells = as.numeric(clust_counts[!!sym(cluster_col)])
    )

    res <- dplyr::relocate(res, .data$n_cells, .before = .data$freq)
  }

  res <- dplyr::mutate(res, pct = (.data$freq / .data$n_cells) * 100)

  res
}


#' Summarize values for chains
#'
#' Summarize values present for each column provided to the data_cols argument.
#' For each cell, the function(s) provided will be applied to each unique label
#' in chain_col.
#'
#' @param sobj_in Seurat object containing V(D)J data
#' @param data_cols meta.data columns to summarize
#' @param fn Function to use for summarizing data_cols
#' @param chain_col meta.data column(s) containing labels for each chain
#' expressed in the cell. These labels are used for grouping the summary
#' output. Set chain_col to NULL to group solely based on the cell ID.
#' @param include_cols Additional columns to include in the output data.frame
#' @param sep Separator to use for expanding data_cols and chain_col
#' @return data.frame containing summary results
#' @export
summarize_chains <- function(sobj_in, data_cols = c("umis", "reads"), fn,
                             chain_col = "chains", include_cols = NULL, sep = ";") {

  # Fetch meta.data
  fetch_cols <- c(data_cols, chain_col, include_cols)

  meta_df <- Seurat::FetchData(sobj_in, fetch_cols)
  meta_df <- tibble::as_tibble(meta_df, rownames = ".cell_id")

  meta_df <- dplyr::filter(
    meta_df,
    across(all_of(data_cols), ~ !is.na(.x))
  )

  # Expand meta.data
  res <- dplyr::mutate(meta_df, across(
    all_of(c(data_cols, chain_col)),
    ~ strsplit(as.character(.x), sep)
  ))

  res <- tidyr::unnest(res, cols = all_of(c(data_cols, chain_col)))

  # Summarize data_cols for each chain present for the cell
  res <- dplyr::mutate(res, across(
    all_of(data_cols),
    ~ .convert_char(.x, as.numeric)
  ))

  grp_cols <- c(".cell_id", chain_col, include_cols)
  res      <- dplyr::group_by(res, !!!syms(grp_cols))

  res <- dplyr::summarize(
    res,
    across(all_of(data_cols), fn),
    .groups = "drop"
  )

  res
}




#' Split V(D)J meta.data columns into vectors
#'
#' @param df_in data.frame to modify
#' @param sep Separator to use for expanding meta.data columns
#' @param sep_cols Columns to split based on sep
#' @param num_cols Columns to convert to numeric
#' @param lgl_cols Columns to convert to logical
#' @return data.frame
.split_vdj <- function(df_in, sep = ";", sep_cols, num_cols = c("umis", "reads"),
                       lgl_cols = c("productive", "full_length")) {

  # Split columns into vectors
  res <- dplyr::mutate(df_in, !!!syms(sep_cols))

  res <- dplyr::mutate(res, across(
    all_of(unname(sep_cols)),
    ~ strsplit(as.character(.x), sep)
  ))

  # Convert to numeric or logical
  res <- dplyr::mutate(
    res,
    across(
      all_of(num_cols),
      map, ~ .convert_char(.x, as.numeric)
    ),
    across(
      all_of(lgl_cols),
      map, ~ .convert_char(.x, as.logical)
    )
  )

  res
}


#' Identify columns with V(D)J data
#'
#' @param df_in data.frame
#' @param clone_col Column containing clonotype IDs to use for identifying
#' columns with V(D)J data
#' @param cols_in meta.data columns containing VDJ data to use for filtering.
#' If set to NULL (recommended) columns are automatically selected.
#' @param sep Separator to search for in columns
#' @return data.frame
.get_vdj_cols <- function(df_in, clone_col, cols_in, sep) {

  # Identify columns with VDJ data based on NAs in clonotype_col
  # If no clonotype_col is given use all columns
  if (is.null(clone_col)) {
    cols_in <- colnames(df_in)
    cols_in <- cols_in[cols_in != ".cell_id"]
  }

  if (is.null(cols_in)) {
    cols_in <- dplyr::mutate(
      df_in,
      across(dplyr::everything(), is.na)
    )

    cols_in <- purrr::keep(
      cols_in,
      ~ identical(.x, pull(cols_in, clone_col))
    )

    cols_in <- colnames(cols_in)
  }

  # Identify columns to split based on sep
  sep_cols <- NULL

  if (!is.null(sep)) {
    sep_cols <- dplyr::select(df_in, all_of(cols_in))

    sep_cols <- purrr::keep(
      sep_cols,
      ~ any(purrr::map_lgl(na.omit(.x), grepl, pattern = sep))
    )

    sep_cols <- colnames(sep_cols)
  }

  # Return list of vectors
  res <- list(
    "vdj" = cols_in,
    "sep" = sep_cols
  )

  res
}


#' Attempt to convert character vector using provided function
#'
#' @param x Character vector to convert
#' @param fn Function to try
#' @return Value converted using fn
.convert_char <- function(x, fn) {
  if (!is.character(x)) {
    return(x)
  }

  suppressWarnings(ifelse(!is.na(fn(x)) | is.na(x), fn(x), x))
}
