## Single-cell Pipeline for Protocol 23_301 (APECED Trial)

The following figure summarizes the computational workflow used to process duodenal biopsy samples from **four paired subjects** and one spike control sample collected **pre- and post-JAK inhibition**. Samples were digested, pooled, sorted, and processed using **10X Chromium GEM-X v4 chemistry**, followed by demultiplexing and downstream single-cell analysis.

### Summary

* **Primary analysis** was performed on paired duodenal biopsies collected before and after ruxolitinib treatment as part of the APECED clinical trial.
* **Cell Ranger** (GEX1–6) generated raw single-cell outputs, followed by doublet detection and demultiplexing.
* **Integration and annotation** used the [Pan-GI Gut Cell Atlas](https://www.gutcellatlas.org/pangi.html) and also Monaco reference annotation .
* A **seurat split** produced lymphoid and non-lymphoid subsets for downstream analyses.
* Final steps included **limma** and **GSEA** (and **CellChat**), with all results compiled into HTML reports.

### Study Aim

Autoimmune PolyEndocrinopathy-Candidiasis-Ectodermal Dystrophy (**APECED**) is caused by loss-of-function mutations in **AIRE**, resulting in multi-organ autoimmunity. Since disease pathology may involve chronic **type 1 inflammatory signaling**, a clinical trial is testing **ruxolitinib** to inhibit the **JAK-STAT pathway**.
By profiling **oral mucosal and intestinal tissues** at single-cell resolution before and after treatment, this study seeks to characterize immune dysregulation in APECED and evaluate cellular responses to JAK inhibition.

## RUN LOCATION
/data/chi/PROJECTS/23-306_Manthiram/Bioinformatics/Runs/Flowcell_22CFHGLT1_Flowcell_232LTWLT4/cellranger_out

### **Pipeline V2 Overview**

<div align="center">
<img src="https://github.niaid.nih.gov/raw/cheungf/23_301/refs/heads/main/IMAGES/pipeline.png?token=GHSAT0AAAAAAAAAA24LDKYKAWG7VWX5W53I2KJRS6A" width="800"/>
</div>


**Figure:** High-level workflow summarizing lineage splitting, atlas integration, spike QC, marker comparison, and downstream UMAP assessments.

---

### **Step 4B - Spike UMAP**

<div align="center">
<img src="https://github.niaid.nih.gov/raw/cheungf/23_301/refs/heads/main/IMAGES/spike_batch.png?token=GHSAT0AAAAAAAAAA24KVOF2YYD3RXZAFL7U2KJRVUA" width="800"/>
</div>


**Figure:** Spike-in samples distribute according to biological structure rather than sequencing batch, indicating successful batch alignment and QC.

---

### **Step 5B - Pan-GI Atlas vs. Single-Gene Markers**

<div align="center">
<img src="https://github.niaid.nih.gov/raw/cheungf/23_301/refs/heads/main/IMAGES/PanGI_L3_violin_markers_marker_celltype.png?token=GHSAT0AAAAAAAAAA24KHYC2Z4EWJL3MSW4M2KJRWMQ" width="800"/>
</div>


**Figure:** Comparison of PanGI level-3 annotations with established single-gene markers shows consistent lineage-level assignments and validates the atlas-driven classification.

---

## **Step 6 - Lineage Split (CD45+ vs CD45–)**

### **Rationale**

The PanGI atlas contains both immune and tissue-derived lineages. Splitting cells by **CD45 (PTPRC)** provides a biologically grounded separation of **Haemopoietic** (immune) versus **Non-Haemopoietic** (tissue/stromal) compartments before downstream analyses.

### **Summary of Findings**

* **Full PanGI UMAP:** CD45 expression cleanly separates immune vs. non-immune clusters across the whole object.
* **CD45+ (Haemopoietic) Subset:** Contains T, B, NK, and myeloid clusters, with clear immune-lineage structure.
* **CD45– (Non-Haemopoietic) Subset:** Contains epithelial, endothelial, stromal, and neural populations as expected.
* **Mitochondrial content:** Highest in the CD45– subset, consistent with epithelial/stromal metabolic activity.
* **Monaco annotations:** Included as an external reference; while coarse, they correctly reflect immune vs. non-immune grouping.

### **Key Takeaway**

**CD45 expression, PanGI lineage annotations, and mitochondrial levels all reinforce a biologically consistent lineage separation, validated across independent marker systems.**

---

### **UMAPs**

---

### **Before Lineage Split (Full PanGI Object)**

<div align="center">
<img src="https://github.niaid.nih.gov/raw/cheungf/23_301/refs/heads/main/IMAGES/BEFORE_SPLIT_UMAP_Annotations_CD45_percentMito.png?token=GHSAT0AAAAAAAAAA24LRZWMHQ42XIWWGTA42KJRXCQ" width="800"/>
</div>

**Figure:** Full dataset with PanGI annotations, CD45 expression, and mitochondrial percentage visualized prior to lineage separation.

---

### **CD45+ (Haemopoietic) Subset**

<div align="center">
<img src="https://github.niaid.nih.gov/raw/cheungf/23_301/refs/heads/main/IMAGES/UMAP_Annotations_CD45_percentMito_Haemopoietic.png?token=GHSAT0AAAAAAAAAA24LWBDOBM623NGCFV342JZZAEA" width="800"/>
</div>

**Figure:** Immune-lineage subset displaying T, B, NK, and myeloid populations with expected CD45 and %mito profiles.

---

### **CD45- (Non-Haemopoietic) Subset**

<div align="center">
<img src="https://github.niaid.nih.gov/raw/cheungf/23_301/refs/heads/main/IMAGES/UMAP_Annotations_CD45_percentMito_NonHaemopoietic.png?token=GHSAT0AAAAAAAAAA24KEERIZGRPXXZMYHGM2JZZAFQ" width="800"/>
</div>

**Figure:** Non-immune subset containing epithelial, endothelial, stromal, and neural lineages with characteristically higher %mito signals.

---

### **Overall Conclusion**

Pipeline V2 integrates spike-in QC, atlas-driven annotations, and CD45-based lineage splitting into a coherent workflow. The agreement between **CD45**, **PanGI lineage labels**, and **mitochondrial metrics** demonstrates that the lineage separation is biologically valid and reproducible across independent markers.

---
