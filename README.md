# GSE40012-Data-Normalization-and-Subtype-Analysis
GSE40012 Data Normalization and Subtype Analysis
# 📊 Microarray Preprocessing Pipeline: Distinguishing Acute Inflammation Subtypes (GSE40012)

This R script automates the essential first steps for analyzing the **GSE40012** microarray dataset. This crucial study investigates the transcriptional differences among patients with **Influenza**, **Bacterial Pneumonia**, **SIRS (Non-Infectious)**, and **Healthy Controls**.

The pipeline performs robust data quality control (QC), normalization, low-expression filtering, and preliminary visualization to ensure the data is prepared for comparative statistical analysis.

## 🚀 Key Features

* **Automated Data Retrieval:** Fetches expression data and metadata directly from the **GEO database (GSE40012)** using `GEOquery`.
* **Quantile Normalization:** Applies **quantile normalization** (`limma::normalizeBetweenArrays`) to correct for technical variation and ensure comparable gene distributions across all samples.
* **Complex Subtyping:** Uses the `dplyr::case_when` function to accurately categorize samples into four distinct groups: **Influenza A**, **Bacterial Pneumonia**, **SIRS (Non-Infectious)**, and **Healthy Control**.
* **QC Visualization:** Generates **boxplots** before and after normalization to visually confirm data distribution consistency.
* **Low-Expression Filtering:** Filters out genes with mean expression in the lowest quartile (25th percentile) to reduce noise.
* **Subtype Clustering:** Performs and visualizes **Principal Component Analysis (PCA)** to check for clear sample separation based on the patient's condition.
* **Data Persistence:** Saves the fully processed and filtered expression matrix, metadata, and subtype information to an `.RData` file.

---

## 🔬 Analysis Overview

| Component | Method / Test | Purpose |
| :--- | :--- | :--- |
| **Dataset** | GSE40012 | Gene expression in blood of patients with different causes of acute inflammation. |
| **Normalization** | Quantile Normalization | Standardizes gene expression distributions across all arrays/samples. |
| **Filtering** | Interquartile Range (IQR) Filtering | Removes non-informative, low-expressed genes for improved statistical power. |
| **Clustering** | PCA | Assesses global sample similarity and confirms separation between the four clinical groups. |

---

## 🛠️ Prerequisites and Setup

### 📦 Packages

The script loads the following essential packages:
* `GEOquery` (For data retrieval)
* `limma` (For normalization)
* `dplyr` (For data manipulation and complex subtype extraction)
* `ggplot2` (For PCA visualization)

### ⚙️ Execution

1.  **Download** the `GSE40012 Data Normalization and Subtype Analysis.R` file.
2.  **Execute** the script in your R environment:
    ```R
    source("GSE40012 Data Normalization and Subtype Analysis.R")
    ```
    *Note: All output files are saved to the current working directory where the script is executed.*

---

## 📁 Output Files (3 Plots + 1 Data File)

| Filename | Type | Description |
| :--- | :--- | :--- |
| `GSE40012_processed_data_step1.RData` | R Binary Data | Contains the final, filtered, and normalized `expression_data`, `metadata`, and `subtype` factor for downstream DGE analysis. |
| `GSE40012_boxplot_before_normalization.png` | QC | Boxplot illustrating raw data distributions. |
| `GSE40012_boxplot_after_normalization.png` | QC | Boxplot confirming uniform data distributions across all samples after quantile normalization. |
| `GSE40012_pca_plot.png` | Clustering | **Principal Component Analysis (PCA)** plot, showing sample clustering colored by clinical subtype. |
