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
# to get smiles
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
  options("java.parameters"=c("-Xmx12G"))
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
                  "org.openscience.cdk.qsar.descriptors.molecular.CPSADescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.CarbonTypesDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.ChiChainDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.ChiClusterDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.ChiPathClusterDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.ChiPathDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.EccentricConnectivityIndexDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.FMFDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.FragmentComplexityDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.GravitationalIndexDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.HBondAcceptorCountDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.HBondDonorCountDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.HybridizationRatioDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.IPMolecularLearningDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.KappaShapeIndicesDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.KierHallSmartsDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.LargestChainDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.LargestPiSystemDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.LengthOverBreadthDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.LongestAliphaticChainDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.MDEDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.MannholdLogPDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.MomentOfInertiaDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.PetitjeanNumberDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.PetitjeanShapeIndexDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.RotatableBondsCountDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.RuleOfFiveDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.VABCDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.VAdjMaDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.WHIMDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.WeightDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.WienerNumbersDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.XLogPDescriptor",
                  "org.openscience.cdk.qsar.descriptors.molecular.ZagrebIndexDescriptor")
  pb <-  progress::progress_bar$new(total = length(smiles))
  dplyr::bind_rows(lapply(smiles, calculate_single_smile, which.desc, pb)) %>%
    tibble::rownames_to_column("smile")
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
    rJava::.jcall("java/lang/System","V","gc")
    gc()
    f
  },
  error = function(e) {
    warning(e)
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
