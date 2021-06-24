#' @export
drug_descriptors_types <- function() {
  rcdk::get.desc.categories()
}

#' @export
drug_descriptors_names <- function(type = "all") {
  # assert type
  rcdk::get.desc.names(type)
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
  f<- f %>% mutate(across(everything(), .fns = ~tidyr::replace_na(.,0)))
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
    return(NA)
  })
}

