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
  if (!all(!missing(category),
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

#' @export
calculate_drugs_AIO <- function(smiles, id_column, smile_column, which.desc) {
  #asserts and dynamic columns
  f <- data.frame()
  sm <- smiles %>% dplyr::select(!!rlang::sym(smile_column)) %>% pull()
  id <- smiles %>% dplyr::select(!!rlang::sym(id_column)) %>% pull()
  for (i in 1:nrow(smiles)) {
    f <- f %>% dplyr::bind_rows(calculate_smile_features(sm[i],
                                                         id[i],
                                                         which.desc))
  }
  f[ifelse(sapply(f, function(x)all(is.na(x))) == TRUE, TRUE, FALSE)] <- NULL
  f[ifelse(sapply(f, function(x)all(is.nan(x))) == TRUE, TRUE, FALSE)] <- NULL
  f[ifelse(sapply(f, function(x)all(x==0)) == TRUE, TRUE, FALSE)] <- NULL
  # f<- f %>% mutate(across(everything(), .fns = ~tidyr::replace_na(.,0)))
  rownames(f) <- f[["drug"]]
  f
}

calculate_smile_features <- function(smile, drug_id, which.desc) {
  tryCatch({
    # check for ChemminerR existence?
    sdf <- ChemmineR::smiles2sdf(smile)
    sdf <- ChemmineR::generate3DCoords(sdf)
    tmpdest <- tempfile(pattern='xxx')
    ChemmineR::write.SDF(sdf, tmpdest)
    farr <- rJava::.jarray(tmpdest, contents.class = 'S')
    molecules <- rJava::.jcall('org/guha/rcdk/util/Misc',
                               '[Lorg/openscience/cdk/interfaces/IAtomContainer;',
                               'loadMolecules',
                               farr,
                               TRUE,
                               TRUE,
                               TRUE,
                               check = FALSE)
    cbind(drug = drug_id, smile, rcdk::eval.desc(molecules, which.desc))
  },

  error = function(e) {
    NA
  })
}

