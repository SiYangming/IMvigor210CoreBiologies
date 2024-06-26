% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{human_gene_signatures}
\alias{human_gene_signatures}
\title{Human gene signatures}
\format{
list with each member representing one gene set.
\describe{
  \item{CD 8 T effector}{CD8 T-effector signature}
  \item{DDR}{DNA Damage Repair signature (based on KEGG)}
  \item{APM}{Antigen Presentation Machinery signature (based on literature)}
  \item{Immune Checkpoint}{Immune checkpoint gene set}
  \item{CC Reg}{Cell Cycle Regulators gene set (based on literature)}
  \item{Fanconi}{Fanconi anemia pathway gene set (based on KEGG)}
  \item{gene19}{pan-tissue fibroblast TGFb response signature (pan-F-TBRS)}
  \item{tcga}{gene set used for TCGA subtyping (based on literature)}
  \item{Cell cycle}{Cell cycle gene set (based on KEGG)}
  \item{Mismatch Repair}{Mismatch Repair gene set (based on KEGG)}
  \item{Homologous recombination}{Homologous recombination gene set (based on KEGG)}
  \item{Nucleotide excision repair}{Nucleotide excision repair gene set (based on KEGG)}
  \item{DNA replication}{DNA replication gene set (based on KEGG)}
  \item{Base excision repair}{Base excision repair gene set (based on KEGG)}
  \item{EMT1}{Epithelial to Mesenchymal Transition signature 1 (based on literature)}
  \item{EMT2}{Epithelial to Mesenchymal Transition signature 2 (based on literature)}
  \item{EMT3}{Epithelial to Mesenchymal Transition signature 3 (based on literature)}
}
}
\source{
Extended Data Table S2
}
\usage{
data(human_gene_signatures)
}
\description{
A dataset containing all gene signatures used in the manuscript
}
\keyword{datasets}
