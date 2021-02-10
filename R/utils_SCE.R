#' Add additional entries to the colData of an SCE
#' @param SCE.in SCE object
#' @param to_add data.frame with additional colData entries
#' @param merge_on Indicate the column of \code{to_add}, which should contain the
#' entries corresponding to the colnames(SCE.in). Default: NULL, which will use
#'  rownames of the \code{to_add} df.
#' @return SCE object with amended colData entries
#' @import SingleCellExperiment
.add_colData <- function(SCE.in, to_add, merge_on = NULL){

    to_add <- as.data.frame(to_add)

    if(is.null(merge_on)){
        cnames <- rownames(to_add)
    }else{
        cnames <- to_add[, merge_on, drop=TRUE]
    }
    if(!all(cnames %in% colnames(SCE.in))){stop("The rownames/merge_on entries
        of the data.frame to be added aren't all part of the SCE")}
    if(!all(colnames(SCE.in) %in% cnames)){message("There are more cells in the
        SCE object than entries in the df; there will be NAs in the resulting colData")}

    cd <- colData(SCE.in)
    cd$MERGE <- colnames(SCE.in)
    to_add$MERGE <- cnames

    cd.out <- merge(cd, to_add, by = "MERGE", all.x = TRUE)
    rownames(cd.out) <- cd.out$MERGE
    cd.out$MERGE <- NULL

    sce.out <- SCE.in
    colData(sce.out) <- DataFrame(cd.out[colnames(sce.out),])
    return(sce.out)

}
