# Interactive QC Filtering for sc/snRNA-seq

This is a Shiny web application for interactive quality control (QC) filtering of single-cell and single-nucleus RNA sequencing (sc/snRNA-seq) data. The app allows users to explore and filter their data based on customizable thresholds for `nCount_RNA`, `nFeature_RNA`, and `percent.mt`.

## 🚀 Features

- Upload data as `.csv` or `.rds` files
- Select dataframes from the R global environment (when running locally)
- Filter by sample, source (raw/filtered), and QC metrics
- Visualize retained and filtered cells with scatter plots
- Summary tables for remaining, filtered, and percent filtered cells

## 🌐 Live App

This app is hosted using shinyapps and runs entirely in your browser (no R server required).

👉 **Access the app here:** https://ahinsu.shinyapps.io/scrnaseq-thresholdselector/ OR https://ahinsu-scrnaseq-thresholdselector.share.connect.posit.cloud



## 🛠️ How to Run Locally

1. Clone the repository:
   ```bash
   git clone https://github.com/ahinsu/scRNAseq-ThresholdSelector.git
   
2. Open RStudio and run Threshold Selector local.R.


