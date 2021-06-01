#' Builds drugs target association matrix from given drugs_targets dataframe
#'
#' The dataframe must contains drugs and targets columns.
#'
#' @param drugs_targets dataframe that contains drugs_targets data
#' @param drugs_column name of drugs column in the given dataset
#' @param  targets_column name of targets column in the given dataset
#' @return dataframe with following sepecifications:
#'          * first column as drugs ids
#'          * targets count other columns (with targets ids as headers)
#'          * values for these columns will be either:
#'            * 1 there is a relation,
#'            * 0 there is no relation.
#' @examples
#' \dontrun{
#'   build_drgus_target_matrix(drugs_targets = dbdataset::Targets_Drug,
#'                             drugs_column = parent_key,
#'                             targets_column = id)
#' }
#' @export
build_drgus_target_matrix <- function(drugs_targets,
                                      drugs_column,
                                      targets_column) {
  assertthat::assert_that(
    !missing(drugs_targets),
    !is.null(drugs_targets),
    "data.frame" %in% class(drugs_targets),
    nrow(drugs_targets) > 0,
    msg = "drugs_targets must be a data.frame that contains"
  )
  assertthat::assert_that(
    !missing(drugs_column),
    !is.null(drugs_column),
    is.character(drugs_column),
    drugs_column != "",
    drugs_column %in% names(drugs_targets),
    msg = paste("drugs_column must be a character value reprsents drugs column",
                "name in drugs_target data.frame")
  )
  assertthat::assert_that(
    !missing(targets_column),
    !is.null(drugs_column),
    is.character(targets_column),
    targets_column != "",
    targets_column %in% names(drugs_targets),
    msg = paste("targets_column must be a character value reprsents drugs column",
                "name in drugs_target data.frame")
  )
}
