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
calculate_drugs_AIO <- function(smiles, which.desc) {
  #asserts and dynamic columns
  f <- data.frame()
  for (i in 1:nrow(smiles)) {
    message(paste0("smile #", i))
    f <- f %>% dplyr::bind_rows(calculate_smile_features(smiles$value[[i]],
                                                          smiles$parent_key[[i]],
                                                          which.desc))
  }
  f
  # f <- dplyr::bind_cols(lapply(which.desc, calculate_single_feature, mol))
  # # drop NA columns (for now)
  # f$geomShape <- NULL
  # # drop 0s only columns
  # zeros <- ifelse(sapply(f, function(x)all(x==0)) == TRUE, T,F)
  # f[zeros] <- NULL
  # # 6. impute NAs with 0 (until I ask Ali about this)
  # f[is.na(f)] <- 0
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

