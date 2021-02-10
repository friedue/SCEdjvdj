#' Import VDJ data to be added to an SCE object
#'
#' @param vdj_dir cellranger VDJ output directories. If a vector of multiple
#' paths is provided, an equal number of cell prefixes must also be provided.
#' If a named vector is given, the names will be used to prefix each cell
#' barcode.
#' @param prefix Prefix to add to new \code{colData} columns
#' @param cell_prefix Prefix to add to cell barcodes
#' @param filter_contigs Only include chains with at least one productive
#' contig
#' @param sep Separator to use for storing per cell clonotype information in
#' the meta.data
#' @param return_SCE TRUE by default. If set to FALSE, only a data.frame will
#' be returned.
#' @return Data.frame with the VDJ info that can be added to the colData of an SCE.
#' @export
import_vdj <- function(vdj_dir, prefix = "", cell_prefix = "",
    filter_contigs = TRUE, sep = ";") {

    # VDJ columns
    count_cols <- c("reads", "umis")

    sep_cols <- c(
        "v_gene",     "d_gene",
        "j_gene",     "c_gene",
        "chains",     "cdr3",
        "cdr3_nt",    count_cols,
        "productive", "full_length"
    )

    vdj_cols <- c(
        "barcode", "clonotype_id",
        sep_cols
    )

    # Check path names
    if (is.null(names(vdj_dir))) {
        if (length(vdj_dir) != length(cell_prefix)) {
            stop("Must provide a cell prefix for each path passed to vdj_dir.")
        }

        names(vdj_dir) <- cell_prefix
    }

    if (any(is.na(names(vdj_dir)))) {
        stop("Cell prefixes must not include NAs.")
    }

    nms <- names(vdj_dir) != "" & !grepl("_$", names(vdj_dir))
    names(vdj_dir)[nms] <- paste0(names(vdj_dir)[nms], "_")

    # Load contigs
    # Check given dir before adding outs to path
    contigs <- "filtered_contig_annotations.csv"

    contigs <- purrr::map_chr(vdj_dir, ~ {
        path <- dplyr::case_when(
            file.exists(file.path(.x, contigs)) ~ file.path(.x, contigs),
            file.exists(file.path(.x, "outs", contigs)) ~ file.path(.x, "outs", contigs)
        )

        if (is.na(path)) {
            stop(paste0(contigs, " not found in ", vdj_dir, "."))
        }

        path
    })

    contigs <- purrr::map(
        contigs,
        readr::read_csv,
        col_types = readr::cols()
    )

    # Add cell prefixes
    contigs <- purrr::imap(contigs, ~ {
        .x <- dplyr::mutate(
            .x,
            barcode      = paste0(.y, .data$barcode),
            clonotype_id = paste0(.y, raw_clonotype_id)
        )

        .x <- dplyr::rename(.x, chains = .data$chain)
    })

    contigs <- dplyr::bind_rows(contigs)

    # Filter for productive contigs
    if (filter_contigs) {
        contigs <- dplyr::filter(contigs, .data$productive, .data$full_length)
    }

    contigs <- dplyr::select(contigs, all_of(vdj_cols))

    # Check if sep is already in sep_cols
    if (any(grepl(sep, contigs[, sep_cols]))) {
        stop(paste0("sep '", sep, "' is already present, select a different seperator."))
    }

    # Sum contig reads and UMIs for chains
    # Some chains are supported by multiple contigs
    grp_cols <- vdj_cols[!vdj_cols %in% count_cols]
    contigs  <- dplyr::group_by(contigs, !!!syms(grp_cols))

    contigs  <- dplyr::summarize(
        contigs,
        across(all_of(count_cols), sum),
        .groups = "drop"
    )

    # Order chains and CDR3 sequences
    # When the rows are collapsed, the cdr3 sequences must be in the same order
    # for every cell. This is required so the cdr3 columns can be used directly
    # as the clonotype ID
    contigs <- dplyr::arrange(
        contigs,
        .data$barcode, .data$chains, .data$cdr3_nt
    )

    # Extract isotypes from c_gene for IGH chain
    iso_pat <- "^IGH[ADEGM]"

    if (any(grepl(iso_pat, contigs$c_gene))) {
        contigs <- dplyr::group_by(contigs, .data$barcode)

        contigs <- dplyr::mutate(
            contigs,
            isotype = list(
                substr(
                    .data$c_gene, 1,
                    attr(regexpr(iso_pat, .data$c_gene), "match.length", exact = TRUE)
                )
            ),
            isotype = map_chr(.data$isotype, ~ {
                isos <- ifelse(all(.x == ""), "", unique(.x[.x != ""]))

                dplyr::case_when(
                    isos == ""        ~ "None",
                    length(isos) > 1  ~ "Multi",
                    length(isos) == 1 ~ isos
                )
            })
        )

        contigs <- dplyr::group_by(contigs, .data$isotype)
    }

    # Collapse chains into a single row for each cell
    # Include isotype and clonotype_id as groups so that they are included in the
    # summarized results
    contigs <- dplyr::group_by(
        contigs,
        .data$barcode, .data$clonotype_id,
        .add = TRUE
    )

    meta_df <- summarize(
        contigs,
        n_chains = n(),
        across(
            all_of(sep_cols),
            ~ paste0(as.character(.x), collapse = sep)
        ),
        .groups = "drop"
    )


    meta_df <- tibble::column_to_rownames(meta_df, "barcode")
    meta_df <- dplyr::rename_with(meta_df, ~ paste0(prefix, .x))

    return(meta_df)
}
