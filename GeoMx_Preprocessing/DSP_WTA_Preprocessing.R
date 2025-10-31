## DSP WTA Preprocessing steps

set.seed(2023)

## Upload required tools
library(NanoStringNCTools)
library(GeomxTools)
library(ggplot2)


## Data files
datadir <- file.path("path-to-files/")

# automatically list files in each directory for use
DCCFiles <- dir(file.path(datadir, "dcc"), pattern = ".dcc$",
                full.names = TRUE, recursive = TRUE)
PKCFiles <- dir(file.path(datadir, "pkc"), pattern = ".pkc$",
                full.names = TRUE, recursive = TRUE)
SampleAnnotationFile <-
  dir(file.path(datadir, "annotation"), pattern = ".xlsx$",
      full.names = TRUE, recursive = TRUE)

# load data
Ovary <-
  readNanoStringGeoMxSet(dccFiles = DCCFiles,
                         pkcFiles = PKCFiles,
                         phenoDataFile = SampleAnnotationFile,
                         phenoDataSheet = "Template",
                         phenoDataDccColName = "Sample_ID",
                         protocolDataColNames = c("aoi", "roi"),
                         experimentDataColNames = c("panel"))

# Shift counts to one
Ovary<- shiftCountsOne(Ovary, useDALogic = TRUE)
# Default QC cutoffs are commented in () adjacent to the respective parameters
# study-specific values were selected after visualizing the QC results in more
# detail below
QC_params <-
  list(minSegmentReads = 1000, # Minimum number of reads (1000)
       percentTrimmed = 80,    # Minimum % of reads trimmed (80%)
       percentStitched = 80,   # Minimum % of reads stitched (80%)
       percentAligned = 75,    # Minimum % of reads aligned (80%)
       percentSaturation = 50, # Minimum sequencing saturation (50%)
       minNegativeCount = 1,   # Minimum negative control counts (10)
       maxNTCCount = 9000,     # Maximum counts observed in NTC well (1000)
       minNuclei = 100,         # Minimum # of nuclei estimated (100)
       minArea = 1000)         # Minimum segment area (5000)
Ovary <-
  setSegmentQCFlags(Ovary, 
                    qcCutoffs = QC_params)        

# Collate QC Results
QCResults <- protocolData(Ovary)[["QCFlags"]]
flag_columns <- colnames(QCResults)
QC_Summary <- data.frame(Pass = colSums(!QCResults[, flag_columns]),
                         Warning = colSums(QCResults[, flag_columns]))
QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
  ifelse(sum(x) == 0L, "PASS", "WARNING")
})
QC_Summary["TOTAL FLAGS", ] <-
  c(sum(QCResults[, "QCStatus"] == "PASS"),
    sum(QCResults[, "QCStatus"] == "WARNING"))
col_by <- "segment"

# Graphical summaries of QC statistics plot function
QC_histogram <- function(assay_data = NULL,
                         annotation = NULL,
                         fill_by = NULL,
                         thr = NULL,
                         scale_trans = NULL) {
  plt <- ggplot(assay_data,
                aes_string(x = paste0("unlist(`", annotation, "`)"),
                           fill = fill_by)) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = thr, lty = "dashed", color = "black") +
    theme_bw() + guides(fill = "none") +
    facet_wrap(as.formula(paste("~", fill_by)), nrow = 4) +
    labs(x = annotation, y = "Segment_tags, #", title = annotation)
  if(!is.null(scale_trans)) {
    plt <- plt +
      scale_x_continuous(trans = scale_trans)
  }
  plt
}

pdf(file="Ovary_DSP_WTA_QC_trimmed.pdf", height= 4, width= 7)
QC_histogram(sData(Ovary), "Trimmed (%)", col_by, 80)
dev.off()
pdf(file="Ovary_DSP_WTA_QC_Stitched.pdf", height= 4, width= 7)
QC_histogram(sData(Ovary), "Stitched (%)", col_by, 80)
dev.off()
pdf(file="Ovary_DSP_WTA_QC_Aligned.pdf", height= 4, width= 7)
QC_histogram(sData(Ovary), "Aligned (%)", col_by, 75)
dev.off()
pdf(file="Ovary_DSP_WTA_QC_Saturated.pdf", height= 4, width= 7)
QC_histogram(sData(Ovary), "Saturated (%)", col_by, 50) +
  labs(title = "Sequencing Saturation (%)",
       x = "Sequencing Saturation (%)")
dev.off()
pdf(file="Ovary_DSP_WTA_QC_area.pdf", height= 4, width= 7)
QC_histogram(sData(Ovary), "area", col_by, 1000, scale_trans = "log10")
dev.off()
pdf(file="Ovary_DSP_WTA_QC_nuclei.pdf", height= 4, width= 7)
QC_histogram(sData(Ovary), "nuclei", col_by, 100)
dev.off()

# Subsetting our dataset has removed samples which did not pass QC
dim(Ovary)

# Generally keep the qcCutoffs parameters unchanged. Set removeLocalOutliers to 
# FALSE if you do not want to remove local outliers
Ovary <- setBioProbeQCFlags(Ovary, 
                            qcCutoffs = list(minProbeRatio = 0.1,
                                             percentFailGrubbs = 20), 
                            removeLocalOutliers = TRUE)

ProbeQCResults <- fData(Ovary)[["QCFlags"]]

# Define QC table for Probe QC
qc_df <- data.frame(Passed = sum(rowSums(ProbeQCResults[, -1]) == 0),
                    Global = sum(ProbeQCResults$GlobalGrubbsOutlier),
                    Local = sum(rowSums(ProbeQCResults[, -2:-1]) > 0
                                & !ProbeQCResults$GlobalGrubbsOutlier))
qc_df

#Subset object to exclude all that did not pass Ratio & Global testing
ProbeQCPassed <- 
  subset(Ovary, 
         fData(Ovary)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
           fData(Ovary)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)
dim(ProbeQCPassed)

Ovary <- ProbeQCPassed 


## Create gene level count data

# Check how many unique targets the object has
length(unique(featureData(Ovary)[["TargetName"]]))


# collapse to targets
target_Ovary <- aggregateCounts(Ovary)
dim(target_Ovary)

## Limit of Quantification

# Define LOQ SD threshold and minimum value
cutoff <- 2
minLOQ <- 2

# Calculate LOQ per module tested
LOQ <- data.frame(row.names = colnames(target_Ovary))
for(module in modules) {
  vars <- paste0(c("NegGeoMean_", "NegGeoSD_"),
                 module)
  if(all(vars[1:2] %in% colnames(pData(target_Ovary)))) {
    LOQ[, module] <-
      pmax(minLOQ,
           pData(target_Ovary)[, vars[1]] * 
             pData(target_Ovary)[, vars[2]] ^ cutoff)
  }
}
pData(target_Ovary)$LOQ <- LOQ


## Filtering

LOQ_Mat <- c()
for(module in modules) {
  ind <- fData(target_Ovary)$Module == module
  Mat_i <- t(esApply(target_Ovary[ind, ], MARGIN = 1,
                     FUN = function(x) {
                       x > LOQ[, module]
                     }))
  LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
}
# ensure ordering since this is stored outside of the geomxSet
LOQ_Mat <- LOQ_Mat[fData(target_Ovary)$TargetName, ]

## Segment Gene Detection

# Save detection rate information to pheno data
pData(target_Ovary)$GenesDetected <- 
  colSums(LOQ_Mat, na.rm = TRUE)
pData(target_Ovary)$GeneDetectionRate <-
  pData(target_Ovary)$GenesDetected / nrow(target_Ovary)

# Determine detection thresholds: 1%, 5%, 10%, 15%, >15%
pData(target_Ovary)$DetectionThreshold <- 
  cut(pData(target_Ovary)$GeneDetectionRate,
      breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
      labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))

# Generally, 5-10% detection is a reasonable segment filtering threshold. However, based on the experimental design (e.g. segment types, size, nuclei) and tissue characteristics (e.g. type, age), these guidelines may require adjustment.

# I change it to 5% for this data
target_Ovary <-
  target_Ovary[, pData(target_Ovary)$GeneDetectionRate >= .05]

dim(target_Ovary)
## Gene Detection Rate
library(scales) # for percent

# Calculate detection rate:
LOQ_Mat <- LOQ_Mat[, colnames(target_Ovary)]
fData(target_Ovary)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
fData(target_Ovary)$DetectionRate <-
  fData(target_Ovary)$DetectedSegments / nrow(pData(target_Ovary))

# Plot detection rate:
plot_detect <- data.frame(Freq = c(1, 5, 10, 20, 30, 50))
plot_detect$Number <-
  unlist(lapply(c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
                function(x) {sum(fData(target_Ovary)$DetectionRate >= x)}))
plot_detect$Rate <- plot_detect$Number / nrow(fData(target_Ovary))
rownames(plot_detect) <- plot_detect$Freq


ggplot(plot_detect, aes(x = as.factor(Freq), y = Rate, fill = Rate)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = formatC(Number, format = "d", big.mark = ",")),
            vjust = 1.6, color = "black", size = 4) +
  scale_fill_gradient2(low = "orange2", mid = "lightblue",
                       high = "dodgerblue3", midpoint = 0.65,
                       limits = c(0,1),
                       labels = scales::percent) +
  theme_bw() +
  scale_y_continuous(labels = scales::percent, limits = c(0,1),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "% of Segments",
       y = "Genes Detected, % of Panel > LOQ")

## RUVg Normalization

library(EDASeq)
library(RUVSeq)
library(DESeq2)

## Get housekeeping genes from nCounter panel
Housekeeping <- c("FCF1",'DHX16','ABCF1','EDC3','ZNF143','HDAC3',
                  'CNOT10','SAP130','AGK','AMMECR1L','SF3A3','COG7',
                  'TMUB2','ZC3H14','DDX50','MRPS5','ZNF346','MTMR14','ERCC3','EIF2B4','TLK2',
                  'NUBP1','USP39','PRPF38A','NOL7','CNOT4')
counts = exprs(target_Ovary)


meta_data <-pData(target_Ovary)

set <- newSeqExpressionSet(as.matrix(round(counts)),
                           phenoData = meta_data)
k=1

set <- RUVg(set, Housekeeping, k=k)

meta_data <-pData(target_Ovary)
pheno <- pData(set[3])
## Add RUVg covariates to limma design
W_1<- pheno$W_1

colnames(meta_data)
Segments <- meta_data$segment
Patients <- meta_data$Patient_ID
TimePoint <- meta_data$TimePoint
Repeats <- meta_data$Repeats
Age<- meta_data$Age_at_diagnosis_med_52
PFS <- meta_data$PFS
Surgical <- meta_data$Surgical_resection_desciption
Tissue <- meta_data$Tissue_site
Pri_met <- meta_data$Pri_met




y <- DGEList(counts =counts, group = Segments)

#include others in the factor model
table(y$samples$group)


y$samples$Segments <- Segments
y$samples$W_1 <- W_1
y$samples$Patients <-Patients
y$samples$TimePoint <- TimePoint 
y$samples$Repeats <- Repeats
y$samples$Age <-Age
y$samples$PFS <- PFS
y$samples$Surgical <-Surgical
y$samples$Tissue <- Tissue
y$samples$Pri_met <-Pri_met

#Filter and convert to logCPM

keep <- filterByExpr(y)
summary(keep)
y <- y[keep, , keep.lib.sizes=FALSE]

#Apply TMM (trimmed mean of M-values) normalization 
# Calculate scaling factors to convert raw library sizes into effective library sizes.
y <- calcNormFactors(y,method = "TMM")

# Create design matrix

design <- model.matrix(~0+W_1+Segments)

colnames(design)

#
fit <- voomLmFit(y, design, block=y$samples$Repeats,  sample.weights = TRUE, plot=TRUE)

normalized_data <- fit$EList$E

## Compare Tumor vs TME
my.contrast <- makeContrasts(SegmentsTumorvsSegmentsTME=SegmentsTumor-SegmentsTME,
                             
                             levels = design)
fit_contrast<- contrasts.fit(fit, my.contrast)
fit_contrast <- eBayes(fit_contrast)
summary(decideTests(fit_contrast))
