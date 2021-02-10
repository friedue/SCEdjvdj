#' Theme for djvdj plotting functions
#'
#' @param ttl_size Size of axis titles
#' @param txt_size Size of axis text
#' @param ln_size Size of axis lines
#' @param txt_col Color of axis text
#' @param ln_col Color of axis lines
#' @return ggplot theme
#' @export
vdj_theme <- function(ttl_size = 12, txt_size = 8, ln_size = 0.5, txt_col = "black",
    ln_col = "grey85") {
    res <- ggplot2::theme(
        strip.background  = ggplot2::element_blank(),
        strip.text        = ggplot2::element_text(size = ttl_size),
        panel.background  = ggplot2::element_blank(),
        legend.background = ggplot2::element_blank(),
        legend.title      = ggplot2::element_text(size = ttl_size),
        legend.key        = ggplot2::element_blank(),
        legend.text       = ggplot2::element_text(size = txt_size, color = txt_col),
        axis.line         = ggplot2::element_line(size = ln_size, color = ln_col),
        axis.ticks        = ggplot2::element_line(size = ln_size, color = ln_col),
        axis.text         = ggplot2::element_text(size = txt_size, color = txt_col),
        axis.title        = ggplot2::element_text(size = ttl_size, color = txt_col)
    )

    res
}


#' Set min and max values for column
#'
#' @param df_in Input data.frame
#' @param ft Name of column containing feature values
#' @param mn Minimum percent rank
#' @param mx Maximum percent rank
#' @return data.frame with modified feature values
.set_lims <- function(df_in, ft, mn = NULL, mx = NULL) {

    if (is.null(mn) && is.null(mx)) {
        return(df_in)
    }

    ft <- sym(ft)

    res <- dplyr::mutate(
        df_in,
        pct = dplyr::percent_rank(!!ft)
    )

    if (!is.null(mn)) {
        res <- dplyr::mutate(
            res,
            !!ft := ifelse(.data$pct > mn, !!ft, NA),
            !!ft := ifelse(.data$pct <= mn, min(!!ft, na.rm = TRUE), !!ft)
        )
    }

    if (!is.null(mx)) {
        res <- dplyr::mutate(
            res,
            !!ft := ifelse(.data$pct < mx, !!ft, NA),
            !!ft := ifelse(.data$pct >= mx, max(!!ft, na.rm = TRUE), !!ft)
        )
    }

    res <- dplyr::select(res, -.data$pct)

    res
}


#' Add regression line to ggplot object
#'
#' @param gg_in ggplot object
#' @param lab_pos Position of correlation coefficient label. Set to NULL to
#' omit label
#' @param lab_size Size of label
#' @param ... Additional arguments to pass to geom_smooth
#' @return ggplot object with added regression line
.add_lm <- function(gg_in, lab_pos = NULL, lab_size = 3.5, ...) {

    # Add regression line
    res <- gg_in +
        ggplot2::geom_smooth(
            method   = "lm",
            formula  = y ~ x,
            se       = F,
            color    = "black",
            size     = 0.5,
            linetype = 2,
            ...
        )

    # Calculate correlation
    if (!is.null(lab_pos)) {
        gg_df <- gg_in$data

        x <- as_name(gg_in$mapping$x)
        y <- as_name(gg_in$mapping$y)

        gg_df <- dplyr::mutate(
            gg_df,
            r       = broom::tidy(stats::cor.test(!!sym(x), !!sym(y)))$estimate,
            r       = round(.data$r, digits = 2),
            pval    = broom::tidy(stats::cor.test(!!sym(x), !!sym(y)))$p.value,
            cor_lab = paste0("r = ", .data$r, ", p = ", format(.data$pval, digits = 2)),
            min_x   = min(!!sym(x)),
            max_x   = max(!!sym(x)),
            min_y   = min(!!sym(y)),
            max_y   = max(!!sym(y)),
            lab_x   = (.data$max_x - .data$min_x) * lab_pos[1] + .data$min_x,
            lab_y   = (.data$max_y - .data$min_y) * lab_pos[2] + .data$min_y
        )
    }

    # Add correlation coefficient label
    res <- res +
        ggplot2::geom_text(
            data          = gg_df,
            mapping       = ggplot2::aes(.data$lab_x, .data$lab_y, label = .data$cor_lab),
            color         = "black",
            size          = lab_size,
            check_overlap = T,
            show.legend   = F
        )

    res
}


#' Add list of aes params to ggplot object
#'
#' @param gg_in ggplot object
#' @param add_aes List of aes params to add or override
#' @param lay Layer number to modify
#' @return ggplot object
.add_aes <- function(gg_in, add_aes, lay) {

    # Need to use colour instead of color
    aes_names      <- names(add_aes)
    names(add_aes) <- replace(aes_names, aes_names == "color", "colour")

    # Add aes params
    curr_aes <- gg_in$layers[[lay]]$aes_params

    curr_aes[names(add_aes)]       <- add_aes
    gg_in$layers[[lay]]$aes_params <- curr_aes

    gg_in
}


#' Create ggplot heatmap
#'
#' @param df_in data.frame
#' @param x Variable to plot on the x-axis
#' @param y Variable to plot on the y-axis
#' @param .fill Variable to use for the fill color
#' @param clrs Vector of colors for plotting
#' @param na_color Color to use for missing values
#' @param ttl Legend title
#' @param ang Angle of x-axis text
#' @param hjst Horizontal justification for x-axis text
#' @param ... Additional arguments to pass to ggplot
#' @return ggplot object
.create_heatmap <- function(df_in, x, y, .fill, clrs = NULL, na_color = "white",
    ttl = .fill, ang = 45, hjst = 1, ...) {

    if (is.null(clrs)) {
        clrs <- "#6A51A3"
    }

    if (length(clrs) == 1) {
        clrs <- c("grey90", clrs)
    }

    if (!is.null(x)) {
        res <- ggplot2::ggplot(
            df_in,
            ggplot2::aes(!!sym(x), !!sym(y), fill = !!sym(.fill))
        )

    } else {
        res <- ggplot2::ggplot(
            df_in,
            ggplot2::aes("sample", !!sym(y), fill = !!sym(.fill))
        )
    }

    res <- res +
        ggplot2::geom_tile(...) +
        ggplot2::guides(fill = ggplot2::guide_colorbar(title = ttl)) +
        ggplot2::scale_fill_gradientn(colors = clrs, na.value = na_color) +
        vdj_theme() +
        ggplot2::theme(
            axis.title  = ggplot2::element_blank(),
            axis.line   = ggplot2::element_blank(),
            axis.ticks  = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(angle = ang, hjust = hjst)
        )

    res
}


#' Create ggplot bar graph
#'
#' @param df_in data.frame
#' @param x Variable to plot on x-axis
#' @param y Variable to plot on y-axis
#' @param .fill Variable to use for fill color
#' @param clrs Vector of colors for plotting
#' @param y_ttl Title for y-axis
#' @param ang Angle of x-axis text
#' @param hjst Horizontal justification for x-axis text
#' @param ... Additional arguments to pass to ggplot
#' @return ggplot object
.create_bars <- function(df_in, x, y, .fill, clrs = NULL, y_ttl = y, ang = 45,
    hjst = 1, ...) {

    # Reverse bar order
    lvls  <- rev(levels(pull(df_in, x)))
    df_in <- .set_lvls(df_in, x, lvls)

    if (!is.null(.fill)) {
        res <- ggplot2::ggplot(
            df_in,
            ggplot2::aes(!!sym(x), !!sym(y), fill = !!sym(.fill))
        )

    } else {
        res <- ggplot2::ggplot(
            df_in,
            ggplot2::aes(!!sym(x), !!sym(y))
        )
    }

    res <- res +
        ggplot2::geom_col(..., position = "dodge") +
        ggplot2::labs(y = y_ttl) +
        vdj_theme() +
        ggplot2::theme(
            axis.title.x = ggplot2::element_blank(),
            axis.text.x  = ggplot2::element_text(angle = ang, hjust = hjst)
        )

    if (!is.null(clrs)) {
        res <- res +
            ggplot2::scale_fill_manual(values = clrs)
    }

    res
}


#' Set column levels
#'
#' @param df_in data.frame
#' @param clmn Column to modify
#' @param lvls Levels
#' @return data.frame
.set_lvls <- function(df_in, clmn, lvls) {

    if (!is.null(lvls) && !is.null(clmn)) {
        dat <- pull(df_in, clmn)

        if (!is.character(dat) && !is.factor(dat)) {
            warning("Plot levels were not modified, levels are only modified for characters and factors.")

            return(df_in)
        }

        if (!all(pull(df_in, clmn) %in% lvls)) {
            stop(paste0("Not all labels in ", clmn, " are included in plot_levels."))
        }

        df_in <- dplyr::mutate(
            df_in,
            !!sym(clmn) := factor(!!sym(clmn), levels = unique(lvls))
        )
    }

    df_in
}

