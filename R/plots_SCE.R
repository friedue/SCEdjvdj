#' Plot clonotype abundance
#'
#' @param SCE_in Seurat object containing V(D)J data
#' @param clonotype_col meta.data column containing clonotype IDs to use for
#' calculating clonotype abundance
#' @param cluster_col meta.data column containing cluster IDs to use for
#' grouping cells when calculating clonotype abundance
#' @param calc_abundances Default: FALSE. Set to TRUE if the clonotype abundances
#' are not part of the SCE_in yet.
#' @param plot_type Type of plot to create, can be 'bar' or 'line'
#' @param yaxis Units to plot on the y-axis, either "frequency" or "percent"
#' @param plot_colors Character vector containing colors for plotting
#' @param plot_lvls Character vector containing levels for ordering
#' @param label_col meta.data column containing labels to use for plot
#' @param n_clonotypes Number of clonotypes to label
#' @param color_col meta.data column to use for coloring bars
#' @param label_aes Named list providing additional label aesthetics (color,
#' size, etc.)
#' @param facet_rows The number of facet rows. Use this argument if plot_type =
#' 'bar'
#' @param facet_scales If plot_type = 'bar', this argument passes a scales
#' specification to facet_wrap, can be "fixed", "free", "free_x", or "free_y"
#' @param ... Additional arguments to pass to geom_line
#' @return ggplot object
#' @export
plot_abundance <- function(SCE_in, clonotype_col = "cdr3_nt", cluster_col = NULL,
    calc_abundances = FALSE,
    plot_type = "bar",
    yaxis = "percent", plot_colors = NULL, plot_lvls = NULL, label_col = "cdr3",
    n_clonotypes = 5, color_col = NULL, label_aes = list(), facet_rows = 1,
    facet_scales = "free_x", ...) {

    if (!yaxis %in% c("frequency", "percent")) {
        stop("yaxis must be either 'frequency' or 'percent'.")
    }

    if (!plot_type %in% c("bar", "line")) {
        stop("plot_type must be either 'bar' or 'line'.")
    }

    if (plot_type == "bar" && is.null(label_col)) {
        stop("Must include label_col when plot_type = 'bar'.")
    }

    if (is.null(color_col)) {
        color_col <- cluster_col
    }

    # Calculate clonotype abundance
    # Return seurat since label_col is needed from meta.data
    if(!all(c("clone_freq", "clone_pct","clone_shared" ) %in% names(colData(SCE_in)))){calc_abundances <- TRUE}
if(calc_abundances == TRUE){
    SCE_in <- calc_abundance(
        SCE_in,
        clonotype_col = clonotype_col,
        cluster_col   = cluster_col,
        prefix        = ".",
        return_SCE = TRUE
    )}

data_col <- "clone_pct"

if (yaxis == "frequency") {
    data_col <- "clone_freq"
}

meta_df <- colData(SCE_in)
meta_df <- tibble::as_tibble(meta_df, rownames = ".cell_id")
meta_df <- dplyr::filter(meta_df, !is.na(!!sym(clonotype_col)))

abund_cols <- c(
    cluster_col, clonotype_col, label_col,
    data_col, color_col
)

meta_df <- dplyr::select(
    meta_df,
    all_of(abund_cols)
)

meta_df <- dplyr::distinct(meta_df)

# Rank by abundance
if (!is.null(cluster_col)) {
    meta_df <- .set_lvls(meta_df, cluster_col, plot_lvls)
    meta_df <- dplyr::group_by(meta_df, !!sym(cluster_col))
}

meta_df <- dplyr::mutate(
    meta_df,
    rank = dplyr::row_number(dplyr::desc(!!sym(data_col)))
)

# Identify top clonotypes
top_genes <- dplyr::slice_min(
    meta_df,
    order_by  = rank,
    n         = n_clonotypes,
    with_ties = FALSE
)

meta_df <- dplyr::ungroup(meta_df)

# Create bar graph
if (plot_type == "bar") {
    top_genes <- dplyr::arrange(top_genes, !!sym(data_col))

    top_genes <- .set_lvls(
        df_in = top_genes,
        clmn  = label_col,
        lvls  = dplyr::pull(top_genes, label_col)
    )

    res <- .create_bars(
        df_in = top_genes,
        x     = label_col,
        y     = data_col,
        y_ttl = yaxis,
        .fill = color_col,
        clrs  = plot_colors,
        ang   = 45,
        hjst  = 1
    )

    if (!is.null(cluster_col)){
        res <- res +
            ggplot2::facet_wrap(
                stats::as.formula(paste0("~ ", cluster_col)),
                nrow   = facet_rows,
                scales = facet_scales
            )
    }

    return(res)
}

# Plot abundance vs rank
res <- ggplot2::ggplot(meta_df, ggplot2::aes(rank, !!sym(data_col))) +
    ggplot2::labs(y = yaxis) +
    vdj_theme()

if (is.null(color_col)) {
    res <- res +
        ggplot2::geom_line(...)

} else {
    res <- res +
        ggplot2::geom_line(ggplot2::aes(color = !!sym(color_col)), ...)
}

if (!is.null(plot_colors)) {
    res <- res +
        ggplot2::scale_color_manual(values = plot_colors)
}

# Add labels
if (!is.null(label_col)) {
    res <- res +
        ggrepel::geom_text_repel(
            ggplot2::aes(label = !!sym(label_col)),
            data          = top_genes,
            nudge_x       = 500,
            direction     = "y",
            segment.size  = 0.2,
            segment.alpha = 0.2,
            size          = 3
        )

    res <- .add_aes(res, label_aes, 2)
}

res
    }
