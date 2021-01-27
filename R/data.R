#' Count data for the transcriptomes of the participants of the IMvigor210 trial
#'
#' A dataset containing the transcriptomes and sample annotation for the majority of
#' participants of the IMvigor210 trial (NCT02108652/NCT02951767).
#'
#' @format A countDataSet raw counts for all genes as well as basic feature and
#' all sample annotations reported in the manuscript:
#' \describe{
#'   \item{ANONPT_ID}{Identifier, mock patient ID}
#'   \item{Best Confirmed Overall Response}{response to treament, by RECISTv1.1 and independent radiology facility; four levels: progressive disease (PD), stable disease (SD), partial response (PR) and complete response (CR)}
#'   \item{binaryResponse}{binarized response, based on 'Best Confirmed Overall Response', grouping PD/SD patients into one group, CR/PR patients into another}
#'   \item{IC Level}{level of IHC-assessed PD-L1 staining on immune cells (IC), with IC0 meaning <1 percent, IC1 >= 1 percent but <5 percent, IC2+ >= 5 percent immune cells staining for PD-L1; sample-level}
#'   \item{os}{overall survival in months}
#'   \item{censOS}{censoring information for surival data, with 0=alive and 1=dead}
#'   \item{Lund}{cluster membership for Lund molecular subtypes}
#'   \item{Lund2}{Lund molecular subtypes}
#'   \item{TCGA Subtype}{TCGA molecular subtypes}
#'   \item{Enrollment IC}{Representative IC level of a patient that was used during trial enrollment; patient-level}
#'   \item{Immune phenotype}{tumor immune phenotype, based on IHC assessment; one of the three: inflamed - robust CD8+ T-cell infiltration and PD-L1 expression, excluded - T cells accumulating in the extracellular matrix–rich stroma, desert - paucity of infiltrating lymphocytes within the tumour or surrounding stroma}
#'   \item{Sex}{sex, either M (male) or F (female)}
#'   \item{Race}{race}
#'   \item{Intravesical BCG administered}{(Y)es or (N)o}
#'   \item{Baseline ECOG Score}{ranging from 0-2; explanation at http://ecog-acrin.org/resources/ecog-performance-status}
#'   \item{Tobacco Use History}{either NEVER, PREVIOUS, or CURRENT}
#'   \item{Met Disease Status}{describing whether metastases affect lymph nodes only (LN Only), affect also visceral tissue (Visceral) or are present in the liver (Liver)}
#'   \item{TC Level}{level of IHC-assessed PD-L1 staining on tumor cells (TC), with TC0 meaning <1 percent, TC1 >= 1 percent but <5 percent, TC2+ >= 5 percent tumor cells staining for PD-L1; sample-level}
#'   \item{FMOne mutation burden per MB}{tumor mutation burden, i.e. number of mutations per MB; as estimated by FMOne panel}
#'   \item{Sample age}{three groups - more than two years: sample collected more than 2 years, 1-2 years: sample collected within 1-2 years before treatment start, (less than) 1 year: sample collected most recently}
#'   \item{Tissue}{tissue origin of sample}
#'   \item{Received platinum}{Has the patient ever received platinum treatment? (Y)es/(N)o}
#'   \item{Sample collected pre-platinum}{If patient has received platinum, was sample collected before treatment? (Y)es/(N)o}
#' }
#' @name cds
#' @source EGAS00001002556
#' @docType data
#' @keywords datasets
#' @usage data(cds)
NULL

#' Targeted mutation data for the participants of the IMvigor210 trial
#'
#' A dataset containing the mutation calls based on FMOne panel and basic sample annotation for
#' participants of the IMvigor210 trial (NCT02108652/NCT02951767).
#'
#' @format An NChannelSet containing sample annotations as well as all known and likely variants,
#' except for chromosomal rearrangements, as reported by Foundation Medicine, Inc.
#' Each variant type is stored in a separate channel.
#' \describe{
#'   \item{amplifications}{copy number alterations of genes that are recurrently amplified in cancer}
#'   \item{deletions}{copy number alterations of genes that are recurrently deleted in cancer}
#'   \item{gains}{non-focal amplifications (> 20MB) recurrently amplified in cancer}
#'   \item{known_short}{short variants (< 49bp long) that are known to have deleterious effects in cancer}
#'   \item{likely_short}{short variants (< 49bp long) that are likely to have deleterious effects in cancer}
#'   \item{ANONPT_ID}{Identifier, mock patient ID}
#'   \item{Best Confirmed Overall Response}{response to treament, by RECISTv1.1 and independent radiology facility; four levels: progressive disease (PD), stable disease (SD), partial response (PR) and complete response (CR)}
#'   \item{binaryResponse}{binarized response, based on 'Best Confirmed Overall Response', grouping PD/SD patients into one group, CR/PR patients into another}
#'   \item{IC Level}{level of IHC-assessed PD-L1 staining on immune cells (IC), with IC0 meaning <1 percent, IC1 >= 1 percent but <5 percent, IC2+ >= 5 percent immune cells staining for PD-L1}
#'   \item{TC Level}{level of IHC-assessed PD-L1 staining on tumor cells (TC), with TC0 meaning <1 percent, TC1 >= 1 percent but <5 percent, TC2+ >= 5 percent tumor cells staining for PD-L1}
#'   \item{FMOne mutation burden per MB}{tumor mutation burden, i.e. number of mutations per MB; as estimated by FMOne panel}
#' }
#' @docType data
#' @keywords datasets
#' @name fmone
#' @usage data(fmone)
NULL

#' Color palettes for reproducing original graphs
#'
#' A dataset containing all color palettes used in original graphs
#' @format list with each member representing one color palette.
#' @docType data
#' @keywords datasets
#' @name color_palettes
#' @usage data(color_palettes)
NULL

#' Human gene signatures
#'
#' A dataset containing all gene signatures used in the manuscript
#'
#' @format list with each member representing one gene set.
#' \describe{
#'   \item{CD 8 T effector}{CD8 T-effector signature}
#'   \item{DDR}{DNA Damage Repair signature (based on KEGG)}
#'   \item{APM}{Antigen Presentation Machinery signature (based on literature)}
#'   \item{Immune Checkpoint}{Immune checkpoint gene set}
#'   \item{CC Reg}{Cell Cycle Regulators gene set (based on literature)}
#'   \item{Fanconi}{Fanconi anemia pathway gene set (based on KEGG)}
#'   \item{gene19}{pan-tissue fibroblast TGFb response signature (pan-F-TBRS)}
#'   \item{tcga}{gene set used for TCGA subtyping (based on literature)}
#'   \item{Cell cycle}{Cell cycle gene set (based on KEGG)}
#'   \item{Mismatch Repair}{Mismatch Repair gene set (based on KEGG)}
#'   \item{Homologous recombination}{Homologous recombination gene set (based on KEGG)}
#'   \item{Nucleotide excision repair}{Nucleotide excision repair gene set (based on KEGG)}
#'   \item{DNA replication}{DNA replication gene set (based on KEGG)}
#'   \item{Base excision repair}{Base excision repair gene set (based on KEGG)}
#'   \item{EMT1}{Epithelial to Mesenchymal Transition signature 1 (based on literature)}
#'   \item{EMT2}{Epithelial to Mesenchymal Transition signature 2 (based on literature)}
#'   \item{EMT3}{Epithelial to Mesenchymal Transition signature 3 (based on literature)}
#' }
#' @source Extended Data Table S2
#' @docType data
#' @keywords datasets
#' @name human_gene_signatures
#' @usage data(human_gene_signatures)
NULL

#' .Rdata file containing mouse tumor region-of-interest (ROI) perimeter coordinates and
#'  CD3+ immune cell coordinates from three one of three mouse studies (666)
#'
#' @format Each .rdata file represents a single mouse study and contains a hyperframe where
#’ each row represents a mouse.
#' The columns of the hyperframe are:
#' \describe{
#'   \item{Barcode}{barcode label for the mouse slide}
#'   \item{Group}{treatment group identifier}
#'   \item{window}{a ‘owin’ object as defined in the ‘spatstat’ R package that contains coordinates of the tumor ROI border}
#'   \item{y}{a point process pattern (ppp) object as defined in the ‘spatstat’ R package that contains coordinates of the CD3+ immune cells}
#' }
#' @docType data
#' @keywords datasets
#' @name dat19
#' @usage data(dat19)
NULL


#'  Mouse tumor region-of-interest (ROI) perimeter coordinates and
#'  CD3+ immune cell coordinates from three one of three mouse studies (1430)
#'
#' @format Each .rdata file represents a single mouse study and contains a hyperframe where
#’ each row represents a mouse
#' The columns of the hyperframe are:
#' \describe{
#'   \item{Barcode}{barcode label for the mouse slide}
#'   \item{Group}{treatment group identifier}
#'   \item{window}{a ‘owin’ object as defined in the ‘spatstat’ R package that contains coordinates of the tumor ROI border}
#'   \item{y}{a point process pattern (ppp) object as defined in the ‘spatstat’ R package that contains coordinates of the CD3+ immune cells}
#' }
#' @name dat57
#' @usage data(dat57)
#' @docType data
#' @keywords datasets
NULL

#'  Mouse tumor region-of-interest (ROI) perimeter coordinates and
#'  CD3+ immune cell coordinates from three one of three mouse studies (1436)
#'
#' @format Each .rdata file represents a single mouse study and contains a hyperframe where
#’ each row represents a mouse. The columns of the hyperframe are:
#' \describe{
#'   \item{Barcode}{barcode label for the mouse slide}
#'   \item{Group}{treatment group identifier}
#'   \item{window}{a ‘owin’ object as defined in the ‘spatstat’ R package that contains coordinates of the tumor ROI border}
#'   \item{y}{a point process pattern (ppp) object as defined in the ‘spatstat’ R package that contains coordinates of the CD3+ immune cells}
#' }
#' @name dat25
#' @usage data(dat25)
#' @docType data
#' @keywords datasets
NULL
