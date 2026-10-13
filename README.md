# Spatial SRC project

Research code for kidney spatial transcriptomics, comparing samples labeled **SRC, LN, and Normal** in the notebooks. The analyses use a kidney reference atlas to infer cell composition, identify spatial niches, and examine expression differences across glomerular, periglomerular, and tubulointerstitial regions.

The repository contains 57 main Jupyter notebooks in Python and R, covering 10x Visium processing, pseudobulk differential expression, cell2location deconvolution, CellCharter spatial clustering, pathway activity, scTM topic modeling, and comparisons with public datasets. It is a collection of research analyses without a unified command-line entry point or a locked software environment. Raw data, manual annotations, and most intermediate results are not included.

## Repository structure

| Directory | Analyses |
| --- | --- |
| [01_amin_ng_2024](01_amin_ng_2024/) | Load and inspect the Amin NG 2024 reference data, using a file named with accession GSE211785 |
| [02_10x_visium](02_10x_visium/) | Visium quality control, sample-level pseudobulk, DESeq2 differential expression, isolated spot removal, integration, and gene selection |
| [03_kidney_atlas](03_kidney_atlas/) | Kidney atlas format conversion, reference expression estimation, cell2location cell abundance inference, marker analysis, and NMF |
| [04_spatial_clustering](04_spatial_clustering/) | CellCharter clustering based on cell composition, niche proportion comparisons, marker analysis, and PROGENy pathway activity |
| [07_scTM](07_scTM/) | scTM/STAMP topic modeling, topic markers, and visualization of topic scores across groups |
| [08_gene_list](08_gene_list/) | Custom gene set scoring and spatial visualization |
| [09_analysis_podocyte_niche](09_analysis_podocyte_niche/) | Interactive region annotation, glomerular and surrounding region identification, regional pseudobulk, differential expression, and cell composition comparisons |
| [10_compare_with_public_data](10_compare_with_public_data/) | Comparisons with the local `Jianbo_Qing_JASN_2026` dataset, including pseudobulk, PCA, and overlap of differentially expressed genes |
| [11_compare_with_GSE294965](11_compare_with_GSE294965/) | GSE294965 pseudobulk, joint comparisons with project data, and differential expression |

The original directory numbering is retained; `05_*` and `06_*` analysis directories are absent. Files in `.ipynb_checkpoints/` are historical checkpoints rather than workflow entry points.

## Sample groups

The main Visium workflow uses the following 12 sample identifiers. Group assignments come from `group_dict` in the notebooks. These identifiers do not establish the number of independent participants.

| Group | Sample identifiers |
| --- | --- |
| SRC | `15837`, `15852`, `16512`, `21022`, `8262` |
| Normal | `B23-03794-5`, `B23-05785-7-1`, `B23-05785-7-2` |
| LN | `X1842`, `X1845`, `X1899`, `X1912` |

Some branches use identifiers with region suffixes, such as `15852_1` and `15852_2`. When substituting data, update `sample_list`, group mappings, and metadata indices together.

## Data preparation and paths

Notebooks generally expect their working directory to be the directory containing the notebook. They access files through `../../data`, `../../Data`, and `../../results`. To retain these relative paths, place the data and results directories in the **parent directory of the repository**, for example:

```text
workspace/
├── Spatial_SRC_project/
│   ├── 02_10x_visium/
│   ├── 03_kidney_atlas/
│   └── ...
├── data/
│   ├── 10x_visium/<sample>/
│   ├── Amin_NG_2024/
│   ├── kidney_atlas_Nature2023/
│   └── GSE294965_processed_data.h5ad
├── Data/
│   ├── 10x_visium/<sample>/
│   ├── score_gene_list/
│   └── LN_rnaseq/Jianbo_Qing_JASN_2026/scedata.rds
├── results/
└── manuscript/figures/main_fig1/upset_intersections/
```

Both `data` and `Data` occur in the code. On case-sensitive file systems, either standardize the notebook paths or provide the corresponding directories; duplicating the data is unnecessary. Output locations follow each notebook's `out_dir`, which may differ from its filename.

Main inputs:

| Input | Usage and requirements |
| --- | --- |
| 10x Visium data | Used by `02_10x_visium` and other spatial analysis branches; each sample directory must contain the expression matrix, `spatial/` images, and coordinate files expected by the reader |
| `GSE211785_Susztak_SC_SN_ATAC_merged_PostSCVI_final.h5ad` | Read by `01_amin_ng_2024/01_load_data.ipynb` |
| `GSE183277_Kidney_Healthy-Injury_Cell_Atlas_snCv3_Seurat_03282022.h5Seurat` | Read by `03_kidney_atlas/01_viz_data.ipynb`; downstream reference models use `patient`, `library`, and cell type annotations |
| Full-resolution tissue images | Required for SpatialData import and region annotation in `09_analysis_podocyte_niche` |
| Manual region annotation CSV files | Downstream notebooks read `results/09_analysis_podocyte_niche/02_identify_podocyte_niche/update/<sample>_1.csv`, match spot indices, and use the `annotation` column |
| Gene set and intersection CSV files | External files in `Data/score_gene_list/` and `manuscript/figures/main_fig1/upset_intersections/` |
| `scedata.rds` | Public dataset Seurat object used in `10_compare_with_public_data` |
| `GSE294965_processed_data.h5ad` | Used in `11_compare_with_GSE294965`, with sample and group metadata for pseudobulk aggregation |

These filenames and paths are taken from the code. The repository does not include automated downloads or a complete data preparation workflow.

## Software environment

Jupyter and Python/R kernels configured for the selected analysis branches are required. Existing kernel names include `zl_src`, `zl_cell2location`, `Python 3`, and `R`. These refer to the original environments; select or register suitable kernels on your own system.

The main dependencies below are identified from code imports. Versions are not pinned.

| Purpose | Dependencies |
| --- | --- |
| Python data processing and visualization | `scanpy`, `squidpy`, `anndata`, `numpy`, `pandas`, `scipy`, `matplotlib`, `seaborn`, `scikit-learn` |
| Pseudobulk and statistics | `adpbulk`, `statsmodels`, `statannotations` |
| Spatial models and pathways | `cell2location`, `cellcharter`, `sctm`, `decoupler` |
| Interactive region annotation | `spatialdata`, `spatialdata-io`, `spatialdata-plot`, `napari-spatialdata`, `geopandas`, `shapely`; some notebooks load `jupyter_black` |
| R analysis | `DESeq2`, `apeglm`, `Seurat`, `SeuratDisk`, `harmony`, `sva`, `Matrix` |
| R plotting and utilities | `ggplot2`, `ggrepel`, `cowplot`, `ComplexHeatmap`, `circlize`, `dplyr`, `tidyr`, `glue`, `ArchR` |

Configure environments by analysis branch, particularly separating cell2location from interactive annotation tools. The cell2location spatial models are configured for 30,000 training epochs, so plan computing resources according to dataset size. Napari annotation requires a graphical session. The repository does not provide a `requirements.txt`, Conda environment file, or `renv.lock`; the dependency list is not a validated set of compatible versions.

Once the environment and data are ready, launch Jupyter from the repository root:

```bash
jupyter lab
```

Open the target notebook, confirm its kernel and working directory, and check `in_dir`, `out_dir`, `data_dir`, `anno_dir`, and sample lists before executing cells in order. Use Python's `os.getcwd()` or R's `getwd()` to inspect the working directory.

## Analysis workflow

### 1. Visium quality control and overall expression analysis

Start with [02_10x_visium/01_filter_data.ipynb](02_10x_visium/01_filter_data.ipynb). It reads each sample, retains spots with at least 500 detected genes, and writes one `.h5ad` file per sample.

- `02_create_bulk` → `03_de_analysis`: aggregate raw counts by sample and use DESeq2 to compare `SRC vs Normal` and `LN vs Normal`.
- `04_remove_isolated_spots` → `05_integrate_data` → `06_integrate_data` / `07_integrate_data`: remove isolated spots, merge data, and perform dimensionality reduction and integration analyses in Python or R.
- `08_sel_genes`: organize gene lists, with additional inputs from external intersection CSV files.

### 2. Kidney reference atlas and cell composition

In `03_kidney_atlas`, first run `01_viz_data` → `02_create_anndata` to export the reference Seurat object and convert it to AnnData. Then select the appropriate branch:

- Hierarchical cell types: `03_run_cell2location` → `04_viz_markers` / `05_find_markers` → `06_run_cell2location` → `07_viz_cell2location`.
- Merged cell types: `08_merge_celltypes` → `09_run_cell2location` → `10_viz_cell2location`. This branch supplies the downstream CellCharter analysis.

Spatial mapping reads the initial Visium quality control outputs. Current model settings include `N_cells_per_location=30` and `detection_alpha=20`. The notebooks export the fifth percentile of the cell abundance posterior, `q05_cell_abundance_w_sf`. `11_run_nmf` explores cell composition patterns, and `12_viz_gene` visualizes gene expression in the reference data.

### 3. Spatial niches

In `04_spatial_clustering`, `01_run_cellcharter` reads the merged cell type branch, converts cell abundances to proportions, aggregates spatial neighborhoods, and searches for 10–15 clusters.

Subsequent notebooks examine cluster cell composition, niche proportions across samples, and marker genes. `06_pathway_activity` infers pathway activity using decoupler's PROGENy network and ULM.

### 4. Glomerular and surrounding regions

The `09_analysis_podocyte_niche` workflow includes manual steps:

1. `01_prepare_data` converts Visium data to SpatialData/Zarr. `04_prepare_data` handles coordinate and image alignment for selected samples.
2. `02_identify_podocyte_niche` / `03_identify_podocyte_niche` use Napari to draw regions and export annotations for spots within them. Check sample names, coordinate systems, layer names, and output filenames when running these notebooks.
3. Prepare the downstream `update/<sample>_1.csv` files and standardize annotation labels, including `Glomerular` and `unassigned`. This curation step and the complete annotations are not included in the repository.
4. `05_create_bulk` → `06_de_analysis` / `07_de_analysis`: analyze glomerular and unassigned regions.
5. `08_identify_Periglomerular` labels spatial neighbors of glomerular spots as `Periglomerular`, preserves the original glomerular labels, and relabels remaining unassigned regions as `Tubulointerstitial`. Then run `09_create_bulk` → `10_de_analysis_Periglomerular` / `11_de_analysis_Tubulointerstitial` / `12_plot_heatmap`.
6. `13_compare_cell_type` combines region annotations with the cell composition used upstream of CellCharter to compare cell type proportions across regions.

### 5. Additional analyses and public dataset comparisons

- `07_scTM`: fit 15 topics, export gene weights and topic proportions, identify topic markers, and compare topic scores.
- `08_gene_list`: calculate custom gene set scores and visualize their spatial distributions.
- `10_compare_with_public_data`: generate Visium pseudobulk for this project, prepare public dataset pseudobulk, and perform joint PCA, differential expression, and overlap analysis of differentially expressed genes.
- `11_compare_with_GSE294965`: run `01` → `02` → `03` / `04` to generate public dataset pseudobulk and compare it with project results. `01_create_pseudo_bulk` writes to a directory named `01_viz_data`, which matches the downstream input path.

## Current reproducibility limitations

- `07_scTM` and `08_gene_list` read preprocessed outputs under `results/06_spahdmap/`. The corresponding upstream code is absent and these inputs must be prepared separately.
- In `07_scTM/01_run_sctm.ipynb`, the statement saving `scTM_15topics.h5ad` is commented out, while downstream topic comparisons read that file. Check the object contents and enable saving, or provide an existing result.
- `08_gene_list/01_viz_gene_score.ipynb` overrides the sample list to include only `X1912`. Adjust this setting before generating results for all samples.
- Some manual annotation notebooks retain sample-specific layer and output names. For example, `03_identify_podocyte_niche` reads `X1845` but still saves files named `21022`; check these names before reuse.
- Notebook numbering reflects the original analysis organization. Follow the actual input/output dependencies across branches. Stored figures and execution outputs do not establish reproducibility in a new environment.

## Main outputs

Most outputs are written to `results/<analysis_directory>/<step>/` under the repository's parent directory:

- `.h5ad`: spatial expression after quality control, reference atlas data, cell abundances, clustering results, and pseudobulk objects.
- `.csv` / `.tsv`: count matrices, sample metadata, differential expression results, cell composition, spatial cluster labels, and topic weights.
- `.zarr`: SpatialData objects for interactive tissue image viewing and region annotation.
- `.png` / `.pdf`: spatial maps, PCA/UMAP plots, heatmaps, volcano plots, and group comparisons. Some figures appear only in notebook outputs.

## License

The code is distributed under the [MIT License](LICENSE). Copyright © 2026 Zhijian Li.
