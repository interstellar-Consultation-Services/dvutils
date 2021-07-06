#' Returns drugs descriptors categories that can be used to calculate drugs
#' features.
#'
#' Each category contains list of descriptors
#'
#' @return character vector of available descriptors categories
#' @examples
#' \dontrun{
#'   drug_descriptors_names(drug_descriptors_categories()[1])
#' }
#' @export
drug_descriptors_categories <- function() {
  rcdk::get.desc.categories()
}

#' Returns drugs descriptors names that can be used to calculate drugs
#' features.
#'
#' @param category character value represents one of
#' \code{drug_descriptors_categories}
#' @return character vector of given descriptor category names
#' @examples
#' \dontrun{
#'   # calculate all available drugs features
#'   calculate_drugs_AIO(drug_descriptors_names())
#' }
#' @export
drug_descriptors_names <- function(category = "all") {
  if (!all(!is.na(category),
           !is.null(category),
           is.character(category),
           all(category %in% c("all",drug_descriptors_categories())))) {
    warning(paste("category must be 'all' or value(s) from",
                  "drug_descriptors_categories() returned values.",
                  "Setting default value 'all'"))
    category <- "all"
  }

  rcdk::get.desc.names(category)
}

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
calculate_drugs_AIO <- function(smiles, id_column, smile_column, descriptors) {
  assertthat::assert_that(
    !missing(smiles),
    !is.null(smiles),
    "data.frame" %in% class(smiles),
    nrow(smiles) > 0,
    ncol(smiles) >= 2,
    msg = "smiles must be a dataframe that contains at least two columns"
  )

  assertthat::assert_that(
    !missing(id_column),
    !is.null(id_column),
    is.character(id_column),
    length(id_column) == 1,
    id_column != "",
    id_column %in% names(smiles),
    msg = "id_column must be a character value reprsents a valid smiles column name"
  )

  assertthat::assert_that(
    !missing(smile_column),
    !is.null(smile_column),
    is.character(smile_column),
    length(smile_column) == 1,
    smile_column != "",
    smile_column %in% names(smiles),
    msg ="smile_column must be a character value reprsents a valid smiles column name"
  )
  descriptors_3d <- c("org.openscience.cdk.qsar.descriptors.molecular.CPSADescriptor",
                      "org.openscience.cdk.qsar.descriptors.molecular.GravitationalIndexDescriptor",
                      "org.openscience.cdk.qsar.descriptors.molecular.IPMolecularLearningDescriptor",
                      "org.openscience.cdk.qsar.descriptors.molecular.LengthOverBreadthDescriptor",
                      "org.openscience.cdk.qsar.descriptors.molecular.MomentOfInertiaDescriptor",
                      "org.openscience.cdk.qsar.descriptors.molecular.VABCDescriptor",
                      "org.openscience.cdk.qsar.descriptors.molecular.WHIMDescriptor",
                      "org.openscience.cdk.qsar.descriptors.molecular.PetitjeanShapeIndexDescriptor")

  calculate_3d <- any(descriptors %in% descriptors_3d)
  if (calculate_3d &&
      !"ChemmineR" %in% rownames(installed.packages())) {
    warning(paste("'ChemmineR' package is needed for calculating",
                  "3D molecules data.",
                  "Disabling descriptors that needs 3D data"))
    calculate_3d <- FALSE
  }

  if (calculate_3d) {
    descriptors_1d <- setdiff(descriptors, descriptors_3d)
    descriptors_3d <-descriptors[descriptors %in% descriptors_3d]
  }

  features <- data.frame()
  sm <- smiles %>% dplyr::select(!!dplyr::sym(smile_column)) %>% dplyr::pull()
  id <- smiles %>% dplyr::select(!!dplyr::sym(id_column)) %>% dplyr::pull()
  if (calculate_3d) {
    for (i in 1:nrow(smiles)) {
      features <- features %>%
        dplyr::bind_rows(calculate_smile_features(sm[i], id[i], descriptors_1d))
    }

    features_3d <- data.frame()
    for (i in 1:nrow(smiles)) {
      features_3d <- features_3d %>%
        dplyr::bind_rows(
          calculate_smile_features(sm[i], id[i], descriptors_3d, TRUE))
    }
    features <- merge(features, features_3d)

  } else {
    for (i in 1:nrow(smiles)) {
      features <- features %>%
        dplyr::bind_rows(calculate_smile_features(sm[i], id[i], descriptors))
    }
  }

  features[ifelse(sapply(features, function(x)all(is.na(x))) == TRUE, TRUE, FALSE)] <- NULL
  features[ifelse(sapply(features, function(x)all(is.nan(x))) == TRUE, TRUE, FALSE)] <- NULL
  features[ifelse(sapply(features, function(x)all(x==0)) == TRUE, TRUE, FALSE)] <- NULL
  features[, -c(1,2)] <- useful::simple.impute.data.frame(features[, -c(1,2)])#features %>% mutate(across(everything(), .fns = ~tidyr::replace_na(.,0)))
  rownames(features) <- features[["drug"]]
  features
}

calculate_smile_features <- function(smile,
                                     drug_id,
                                     descriptor,
                                     require_3d = FALSE) {
  tryCatch({
    molecule <- R.utils::withTimeout(get_molecules(smile, require_3d),
                                     timeout = 180)
    cbind(drug = drug_id, smile, rcdk::eval.desc(molecule, descriptor))
  },

  error = function(e) {
    NA
  })
}

get_molecules <- function(smile, require_3d = FALSE) {
  if (require_3d) {
    sdf <- ChemmineR::smiles2sdf(smile)
    sdf <- R.utils::withSink(ChemmineR::generate3DCoords(sdf),
                             file = tempfile())
    tmpdest <- tempfile(pattern='xxx')
    ChemmineR::write.SDF(sdf, tmpdest)
    farr <- rJava::.jarray(tmpdest, contents.class = 'S')
    rJava::.jcall('org/guha/rcdk/util/Misc',
                  '[Lorg/openscience/cdk/interfaces/IAtomContainer;',
                  'loadMolecules',
                  farr,
                  TRUE,
                  TRUE,
                  TRUE,
                  check = FALSE)
  } else {
    rcdk::parse.smiles(smile)
  }
}

