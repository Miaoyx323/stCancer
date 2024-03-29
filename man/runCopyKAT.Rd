% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CNV.R
\name{runCopyKAT}
\alias{runCopyKAT}
\title{runCopyKAT}
\usage{
runCopyKAT(
  SeuratObject,
  savePath,
  cluster.info = "seurat_clusters",
  id.type = "S",
  ngene.chr = 5,
  win.size = 25,
  KS.cut = 0.1,
  n.cores = 1,
  sam.name = "test",
  distance = "euclidean",
  norm.cell.names = "",
  output.seg = "FLASE",
  ...
)
}
\arguments{
\item{SeuratObject}{A seuratObejct to analyze CNV}

\item{savePath}{A path to store the results}

\item{cluster.info}{The label to store cluster information in `SeuratObject@meta.data`}

\item{id.type}{The parameter of copyKAT}

\item{ngene.chr}{The parameter of copyKAT}

\item{win.size}{The parameter of copyKAT}

\item{KS.cut}{The parameter of copyKAT}

\item{n.cores}{The parameter of copyKAT}

\item{sam.name}{The parameter of copyKAT}

\item{distance}{The parameter of copyKAT}

\item{norm.cell.names}{The parameter of copyKAT}

\item{output.set}{The parameter of copyKAT}
}
\value{
A list containing a seuratObject, a list of results of intermediate process,
a list of plots and a logical value to show whether the process completed successfully
}
\description{
Analyze CNV using R package copyKAT
}
