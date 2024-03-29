% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/downsample.R
\name{downsample}
\alias{downsample}
\title{Downsample sce object}
\usage{
downsample(sce, n, seed = 1234)
}
\arguments{
\item{sce}{a SingleCellExperiment object}

\item{n}{integer, number of cells to be downsampled}

\item{seed}{numeric vector; the seed set for each run of clustering to keep the results consistent between runs}
}
\value{
a sce object
}
\description{
Decrease the number of cells per sample by random selection to help reduce the computational burden
which may incur during data processing or clustering
}
\examples{
sce_d <- downsample(sce, n=2000)

}
