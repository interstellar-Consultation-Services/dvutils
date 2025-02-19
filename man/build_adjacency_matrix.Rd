% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools.R
\name{build_adjacency_matrix}
\alias{build_adjacency_matrix}
\title{Builds adjacency matrix from given data between source and target columns}
\usage{
build_adjacency_matrix(data, source, target)
}
\arguments{
\item{data}{dataframe that contains data data}

\item{source}{name of drugs column in the given dataset}

\item{targets}{name of targets column in the given dataset}
}
\value{
dataframe with following sepecifications:
* first column as drugs ids
* targets count other columns (with targets ids as headers)
* values for these columns will be either:
* 1 there is a relation,
* 0 there is no relation.
}
\description{
The dataframe must contains drugs and targets columns.
}
\examples{
\dontrun{
  build_adjacency_matrix(data = dbdataset::Targets_Drug,
                         source = "parent_key",
                         targets = "id")
}
}
