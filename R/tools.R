#' Builds adjacency matrix from given data between source and target columns
#'
#' The dataframe must contains drugs and targets columns.
#'
#' @param data dataframe that contains data data
#' @param source name of drugs column in the given dataset
#' @param targets name of targets column in the given dataset
#' @return dataframe with following sepecifications:
#'          * first column as drugs ids
#'          * targets count other columns (with targets ids as headers)
#'          * values for these columns will be either:
#'            * 1 there is a relation,
#'            * 0 there is no relation.
#' @examples
#' \dontrun{
#'   build_adjacency_matrix(data = dbdataset::Targets_Drug,
#'                          source = "parent_key",
#'                          targets = "id")
#' }
#' @export
build_adjacency_matrix <- function(data,
                                   source,
                                   target) {
  assertthat::assert_that(
    !missing(data),
    !is.null(data),
    "data.frame" %in% class(data),
    nrow(data) > 0,
    msg = "data must be a data.frame that contains"
  )
  assertthat::assert_that(
    !missing(source),
    !is.null(source),
    is.character(source),
    source != "",
    source %in% names(data),
    msg = paste("source must be a character value reprsents drugs column",
                "name in drugs_target data.frame")
  )
  assertthat::assert_that(
    !missing(target),
    !is.null(source),
    is.character(target),
    target != "",
    target %in% names(data),
    msg = paste("target must be a character value reprsents drugs column",
                "name in drugs_target data.frame")
  )
  # not all drugs has targets, what we will do with other drugs?!!!
  # investigate more later
  data %>%
    group_by(!!sym(source), !!sym(target)) %>%
    tally() %>%
    spread(!!sym(target), n, drop = FALSE, fill = 0)
}
