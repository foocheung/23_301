PBMC UMAP QC — Results Summary
================================

WHAT THESE UMAPs SHOW
----------------------
As requested in yesterdays meeting I have generate the required UMAPs:

Three UMAP reductions are compared to separate technical from biological structure:

  - No Harmony         : raw PCA embedding, no batch correction
  - Harmony (pool)     : corrects for sequencing pool only (technical)
  - Harmony (subj+pool): corrects for pool AND inter-subject variation

Each reduction is coloured by Condition, Pool, Monaco cell type, and Subject,
and shown for all samples together as well as each condition group individually
(HC only, Severe only, Treated only, Untreated only).

Per-subject and per-subject x pool subset panels are also included so each
subject's cells are shown in full without background dimming.


KEY FINDINGS
-------------
- Comparing No Harmony vs Harmony (Pool), the UMAPs look very similar, showing
  good mixing of samples across clusters, which indicates little to no strong
  pool-driven batch effect.

- Most cell types are well mixed across samples, but monocytes remain partially
  separated, suggesting this pattern is not technical but reflects real
  biological heterogeneity.

- When looking at control-only and severe-only UMAPs, the same structure is
  preserved, confirming that clustering is not driven by batch but by biological
  differences between subjects.

- The strongest heterogeneity is seen in severe subjects, especially within
  monocytes, where each subject forms slightly distinct subclusters.

- Applying Harmony with Pool + Subject removes this separation in monocytes,
  indicating that this step is correcting subject-specific variation and may
  mask true biological differences.


Side Note:
  DEG and GSEA results are derived entirely from raw count pseudobulk
  models and limma-derived statistics, which are run independently of
  any dimensionality reduction or batch correction applied here.
  Harmony coordinates are never used as input to those models.


