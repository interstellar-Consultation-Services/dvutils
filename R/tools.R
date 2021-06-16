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
    msg = paste("source must be a character value reprsents a valid column",
                "name in data")
  )
  assertthat::assert_that(
    !missing(target),
    !is.null(source),
    is.character(target),
    target != "",
    target %in% names(data),
    msg = paste("target must be a character value reprsents a valid",
                "name in data")
  )
  # not all drugs has targets, what we will do with other drugs?!!!
  # investigate more later
  data %>%
    group_by(!!sym(source), !!sym(target)) %>%
    tally() %>%
    spread(!!sym(target), n, drop = FALSE, fill = 0)
}
# to get smiles
# drugs_with_targets <- unique(dbdataset::Targets_Drug$parent_key)
# smiles <- dbdataset::Calculated_Properties_Drug %>%
#               filter(kind == "SMILES", parent_key %in% drugs_with_targets)
# writeLines(smiles$value[[1]], "test.smi")
# smi = "test.smi"
# mol <- Rcpi::extractDrugAIO(Rcpi::readMolFromSmi(smi))
# there are 440 bioteck drug with smiles
# small_mol <- dbdataset::Drugs %>% filter(type == "small molecule")
# d <- setdiff(drugs_with_targets, small_mol$primary_key)
#' @export
calculate_drugs_AIO <- function(smiles) {
  #options("java.parameters"=c("-Xmx12G"))
  which.desc <- c("org.openscience.cdk.qsar.descriptors.molecular.ALOGPDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.APolDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.AminoAcidCountDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.AromaticAtomsCountDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.AromaticBondsCountDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.AtomCountDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.AutocorrelationDescriptorCharge",
                  "org.openscience.cdk.qsar.descriptors.molecular.AutocorrelationDescriptorMass",
                  "org.openscience.cdk.qsar.descriptors.molecular.AutocorrelationDescriptorPolarizability",
                  "org.openscience.cdk.qsar.descriptors.molecular.BCUTDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.BPolDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.BondCountDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.CarbonTypesDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.ChiChainDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.ChiClusterDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.ChiPathClusterDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.ChiPathDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.EccentricConnectivityIndexDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.FMFDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.FragmentComplexityDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.HBondAcceptorCountDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.HBondDonorCountDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.HybridizationRatioDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.KappaShapeIndicesDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.KierHallSmartsDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.LargestChainDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.LargestPiSystemDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.LongestAliphaticChainDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.MDEDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.MannholdLogPDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.PetitjeanNumberDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.PetitjeanShapeIndexDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.RotatableBondsCountDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.RuleOfFiveDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.VAdjMaDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.WeightDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.WienerNumbersDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.XLogPDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.ZagrebIndexDescriptor")
  # pb <-  progress::progress_bar$new(total = length(smiles))
  # dplyr::bind_rows(lapply(smiles, calculate_single_smile, which.desc, pb)) %>%
  #   tibble::rownames_to_column("smile")
  f <- dplyr::bind_cols(lapply(which.desc, calculate_single_feature, mol))
}

calculate_single_feature <- function(which.desc, mol) {
  message(which.desc)
  rcdk::eval.desc(mol, which.desc)
}

calculate_single_smile <- function(smile, which.desc, pb) {
  pb$tick()
  tryCatch({
    f <- NULL
    if (nchar(smile) < 800) {
      which.desc <- c(which.desc,
                      "org.openscience.cdk.qsar.descriptors.molecular.WeightedPathDescriptor")
    }
    smile <- rcdk::parse.smiles(smile)
    if (!is.null(smile)) {
      f <- rcdk::eval.desc(smile, which.desc, verbose = T)
    }
    if (is.null(f)) {
      f <- data.frame()
    }
    # rJava::.jcall("java/lang/System","V","gc")
    # gc()
    f
  },
  error = function(e) {
    warning(e)
    data.frame()
  })
}

#' @export
calculate_targets_features <- function(targets, id_column, sequence_column) {
  assertthat::assert_that(
    !missing(targets),
    !is.null(targets),
    "data.frame" %in% class(targets),
    nrow(targets) > 0,
    msg = "targets must be a data.frame that contains"
  )
  assertthat::assert_that(
    !missing(id_column),
    !is.null(id_column),
    is.character(id_column),
    id_column != "",
    id_column %in% names(targets),
    msg = paste("id_column must be a character value reprsents targets id column",
                "name in targets dataframe")
  )
  assertthat::assert_that(
    !missing(sequence_column),
    !is.null(sequence_column),
    is.character(sequence_column),
    sequence_column != "",
    sequence_column %in% names(targets),
    msg = paste("sequence_column must be a character value reprsents sequences column",
                "name in targets dataframe")
  )
  # targets <- dbdataset::Polypeptide_Target_Drug %>% dplyr::select(target_id = parent_id, amino_acid_sequence)
  # sequence_column <- "amino_acid_sequence"
  # targets_features <- targets %>% inner_join(features, by = c("amino_acid_sequence", "sequence"))
  new_targets <- targets %>%
    dplyr::select(!!sym(id_column), !!sym(sequence_column))
  new_targets[[sequence_column]] <- sapply(new_targets[[sequence_column]], clean_sequence)
  pb <-  progress::progress_bar$new(total = nrow(new_targets))
  features <- dplyr::bind_rows(apply(new_targets,
                                     1,
                                     calculate_sequence_features,
                                     pb))
}

clean_sequence <- function (sequence) {
  lines <- strsplit(sequence, "\n")[[1]]
  ind <- which(substr(lines, 1L, 1L) == ">")
  nseq <- length(ind)
  if (nseq == 0) {
    stop("no line starting with a > character found")
  }
  start <- ind + 1
  end <- ind - 1
  end <- c(end[-1], length(lines))
  lapply(seq_len(nseq), function(i) paste(lines[start[i]:end[i]],
                                          collapse = ""))[[1]]
}

calculate_sequence_features <- function(target_row, pb) {
    pb$tick()
    sequence <- target_row[2]
    if (!protr::protcheck(sequence)) {
      message(paste("This sequence:",
                    sequence,
                    "is invalid. Skipping Process."))
      return(data.frame())
    }
    tryCatch({
      f <- data.frame(
        taret_id = target_row[1],
        sequence,
        t(protr::extractAAC(sequence)),
        t(protr::extractDC(sequence)),
        t(protr::extractTC(sequence)),
        t(protr::extractCTDC(sequence)),
        #t(protr::extractCTDCClass(sequence)),
        t(protr::extractCTDT(sequence)),
        #t(protr::extractCTDTClass(sequence)),
        t(protr::extractCTDD(sequence)),
        #t(protr::extractCTDDClass(sequence)),
        t(protr::extractCTriad(sequence)),
        #t(protr::extractCTriadClass(sequence)),
        t(protr::extractSOCN(sequence)),
        t(protr::extractQSO(sequence)),
        t(protr::extractPAAC(sequence)),
        t(protr::extractAPAAC(sequence))
      )
      if (length(sequence) > 30) {
        f <- cbind(f, t(protr::extractMoreauBroto(sequence)))
        f <- cbind(f, t(protr::extractMoran(sequence)))
        f <- cbind(f, t(protr::extractGeary(sequence)))
      }
      f
    },
    error = function(e) {
      message(paste("Could not process sequence:",
                    target_row[2],
                    "due to:",
                    e$message))
      data.frame()
    })
  }

# library(dplyr)
# drugs_with_targets <- unique(dbdataset::Targets_Drug$parent_key)
# smiles <- dbdataset::Calculated_Properties_Drug %>%
#   filter(kind == "SMILES", parent_key %in% drugs_with_targets)
#
# mol <- rcdk::parse.smiles(smiles$value[3992], omit.nulls = T)
# f <- rcdk::eval.desc(mol, which.desc, verbose = T)
#
# drugs_features <- smiles %>%
#   select(drug_id = parent_key, smile = value) %>%
#   inner_join(f)

# amin <- dbdataset::Polypeptide_Target_Drug$amino_acid_sequence[[1]]
# protr::extractAAC(protr::readFASTA("test.fasta")[[1]])
