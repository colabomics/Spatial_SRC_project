# Spatial SRC project

用于肾脏空间转录组分析的研究代码，比较代码中标记为 **SRC、LN 和 Normal** 的样本，结合肾脏参考图谱推断细胞组成、识别空间生态位，并分析肾小球、肾小球周围和肾小管间质区域的表达差异。

仓库包含 57 个主要 Jupyter Notebook，使用 Python 和 R，覆盖 10x Visium 数据处理、伪 bulk 差异表达、cell2location 解卷积、CellCharter 空间聚类、通路活性、scTM 主题建模及公共数据比较。当前仓库是研究分析记录，没有统一命令行入口或锁定的运行环境；原始数据、人工标注及大部分中间结果未包含在仓库中。

## 目录导航

| 目录 | 分析内容 |
| --- | --- |
| [01_amin_ng_2024](01_amin_ng_2024/) | 读取并检查 Amin NG 2024 参考数据（文件名包含 GSE211785） |
| [02_10x_visium](02_10x_visium/) | Visium 质控、样本级伪 bulk、DESeq2 差异表达、孤立 spot 移除、数据整合和基因筛选 |
| [03_kidney_atlas](03_kidney_atlas/) | 肾脏图谱格式转换、参考表达谱估计、cell2location 细胞丰度推断、标记基因和 NMF 分析 |
| [04_spatial_clustering](04_spatial_clustering/) | 基于细胞组成的 CellCharter 聚类、生态位比例比较、标记基因和 PROGENy 通路活性 |
| [07_scTM](07_scTM/) | scTM/STAMP 主题建模、主题标记基因和组间主题分数可视化 |
| [08_gene_list](08_gene_list/) | 自定义基因集评分及空间可视化 |
| [09_analysis_podocyte_niche](09_analysis_podocyte_niche/) | 交互式区域标注、肾小球及周围区域识别、区域伪 bulk、差异表达和细胞组成比较 |
| [10_compare_with_public_data](10_compare_with_public_data/) | 与本地 `Jianbo_Qing_JASN_2026` 数据比较，包括伪 bulk、PCA 和差异基因重叠 |
| [11_compare_with_GSE294965](11_compare_with_GSE294965/) | GSE294965 伪 bulk，与本项目数据联合比较及差异表达分析 |

目录保留原始编号；当前仓库没有 `05_*` 和 `06_*` 分析目录。`.ipynb_checkpoints/` 中的文件是历史检查点，不作为运行入口。

## 样本分组

主 Visium 流程使用以下 12 个样本标识，分组来自 Notebook 中的 `group_dict`。这些标识不等同于独立受试者数量。

| 分组 | 样本标识 |
| --- | --- |
| SRC | `15837`、`15852`、`16512`、`21022`、`8262` |
| Normal | `B23-03794-5`、`B23-05785-7-1`、`B23-05785-7-2` |
| LN | `X1842`、`X1845`、`X1899`、`X1912` |

部分分支使用带区域后缀的标识，例如 `15852_1`、`15852_2`；替换数据时需同步检查 `sample_list`、分组映射和元数据索引。

## 数据准备与路径

Notebook 通常以其所在目录为工作目录，通过 `../../data`、`../../Data` 和 `../../results` 访问文件。因此，保留现有相对路径时，数据与结果应位于**仓库的父目录**，例如：

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

`data` 与 `Data` 是代码中实际使用的两种拼写。在区分大小写的文件系统上，需要统一修改 Notebook 的路径，或分别准备对应目录，无需为此重复复制数据。输出位置以各 Notebook 的 `out_dir` 为准，未必与 Notebook 文件名一致。

主要输入包括：

| 输入 | 使用位置与要求 |
| --- | --- |
| 10x Visium 数据 | `02_10x_visium` 和其他空间分析分支；每个样本目录需包含读取函数所需的表达矩阵、`spatial/` 图像及坐标文件 |
| `GSE211785_Susztak_SC_SN_ATAC_merged_PostSCVI_final.h5ad` | `01_amin_ng_2024/01_load_data.ipynb` |
| `GSE183277_Kidney_Healthy-Injury_Cell_Atlas_snCv3_Seurat_03282022.h5Seurat` | `03_kidney_atlas/01_viz_data.ipynb`；下游参考模型使用 `patient`、`library` 和细胞类型注释 |
| 组织全分辨率图像 | `09_analysis_podocyte_niche` 中 SpatialData 导入和区域标注所需 |
| 人工区域标注 CSV | 下游读取 `results/09_analysis_podocyte_niche/02_identify_podocyte_niche/update/<sample>_1.csv`，按 spot 索引匹配并使用 `annotation` 列 |
| 基因集与交集 CSV | `Data/score_gene_list/` 和 `manuscript/figures/main_fig1/upset_intersections/` 中的外部文件 |
| `scedata.rds` | `10_compare_with_public_data` 中公共数据的 Seurat 对象 |
| `GSE294965_processed_data.h5ad` | `11_compare_with_GSE294965`；使用样本及分组元数据构建伪 bulk |

以上文件名和目录来自代码，仓库未提供自动下载或完整数据准备脚本。

## 运行环境

需要 Jupyter，以及按分析分支配置的 Python 和 R 内核。现有内核名称包括 `zl_src`、`zl_cell2location`、`Python 3` 和 `R`，它们只是原运行环境的名称，需在本机重新选择或注册对应内核。

主要依赖根据代码导入整理如下，具体版本未锁定：

| 用途 | 依赖 |
| --- | --- |
| Python 数据处理与可视化 | `scanpy`、`squidpy`、`anndata`、`numpy`、`pandas`、`scipy`、`matplotlib`、`seaborn`、`scikit-learn` |
| 伪 bulk 与统计 | `adpbulk`、`statsmodels`、`statannotations` |
| 空间模型与通路 | `cell2location`、`cellcharter`、`sctm`、`decoupler` |
| 交互式区域标注 | `spatialdata`、`spatialdata-io`、`spatialdata-plot`、`napari-spatialdata`、`geopandas`、`shapely`；部分 Notebook 加载 `jupyter_black` |
| R 分析 | `DESeq2`、`apeglm`、`Seurat`、`SeuratDisk`、`harmony`、`sva`、`Matrix` |
| R 绘图及辅助工具 | `ggplot2`、`ggrepel`、`cowplot`、`ComplexHeatmap`、`circlize`、`dplyr`、`tidyr`、`glue`、`ArchR` |

建议按分支配置环境，尤其将 cell2location 与交互式标注环境分别管理。cell2location 空间模型配置了 30,000 个训练 epoch，应按数据规模准备计算资源；Napari 标注需要图形界面。仓库尚未提供 `requirements.txt`、Conda 环境文件或 `renv.lock`，上述列表不代表已经验证的兼容版本组合。

准备好环境和数据后，从仓库根目录启动：

```bash
jupyter lab
```

打开目标 Notebook，确认内核和工作目录，检查 `in_dir`、`out_dir`、`data_dir`、`anno_dir` 及样本列表，再按单元格顺序执行。可以用 Python 的 `os.getcwd()` 或 R 的 `getwd()` 检查工作目录。

## 分析顺序

### 1. Visium 质控与整体表达分析

从 [02_10x_visium/01_filter_data.ipynb](02_10x_visium/01_filter_data.ipynb) 开始：读取各样本，保留至少检测到 500 个基因的 spot，输出每个样本的 `.h5ad`。

- `02_create_bulk` → `03_de_analysis`：按样本聚合原始计数，使用 DESeq2 比较 `SRC vs Normal` 和 `LN vs Normal`。
- `04_remove_isolated_spots` → `05_integrate_data` → `06_integrate_data` / `07_integrate_data`：移除孤立 spot，合并数据，进行 Python 或 R 中的降维及整合分析。
- `08_sel_genes`：整理基因列表，还依赖仓库外的交集 CSV。

### 2. 肾脏参考图谱与细胞组成

在 `03_kidney_atlas` 中，先运行 `01_viz_data` → `02_create_anndata`，将参考 Seurat 对象导出并转换为 AnnData，然后选择对应分支：

- 分层细胞类型：`03_run_cell2location` → `04_viz_markers` / `05_find_markers` → `06_run_cell2location` → `07_viz_cell2location`。
- 合并细胞类型：`08_merge_celltypes` → `09_run_cell2location` → `10_viz_cell2location`。这是后续 CellCharter 分析读取的分支。

空间映射读取 Visium 初步质控结果。当前模型设置 `N_cells_per_location=30`、`detection_alpha=20`，并导出细胞丰度后验的 5% 分位数 `q05_cell_abundance_w_sf`。`11_run_nmf` 探索细胞组成模式，`12_viz_gene` 展示参考数据中的基因表达。

### 3. 空间生态位

在 `04_spatial_clustering` 中，`01_run_cellcharter` 读取合并细胞类型分支的结果，将细胞丰度转换为比例后进行空间邻域聚合和聚类，搜索 10–15 个簇。

后续 Notebook 分析各簇细胞组成、样本间生态位比例和标记基因；`06_pathway_activity` 使用 decoupler 的 PROGENy 网络与 ULM 推断通路活性。

### 4. 肾小球及周围区域

`09_analysis_podocyte_niche` 的流程包含人工操作：

1. `01_prepare_data` 将 Visium 数据转换为 SpatialData/Zarr；`04_prepare_data` 处理指定样本的坐标及图像对齐。
2. `02_identify_podocyte_niche` / `03_identify_podocyte_niche` 使用 Napari 绘制区域，并将区域内 spot 导出为标注。运行时需核对样本名、坐标系、图层名和输出文件名。
3. 整理下游所需的 `update/<sample>_1.csv`，统一 `Glomerular` 和 `unassigned` 等标注名称。该整理过程及完整标注未包含在仓库中。
4. `05_create_bulk` → `06_de_analysis` / `07_de_analysis`：分析肾小球与未分配区域。
5. `08_identify_Periglomerular` 将肾小球 spot 的空间邻居标为 `Periglomerular`，保留原肾小球标注，并将其余未分配区域标为 `Tubulointerstitial`；随后运行 `09_create_bulk` → `10_de_analysis_Periglomerular` / `11_de_analysis_Tubulointerstitial` / `12_plot_heatmap`。
6. `13_compare_cell_type` 结合区域标注与 CellCharter 上游细胞组成，比较各区域的细胞类型比例。

### 5. 扩展分析与公共数据比较

- `07_scTM`：拟合 15 个主题并导出基因权重、主题比例，提取标记基因并比较主题分数。
- `08_gene_list`：计算自定义基因集分数并绘制空间分布。
- `10_compare_with_public_data`：先生成本项目 Visium 伪 bulk，再准备公共数据伪 bulk，进行联合 PCA、差异表达和差异基因重叠分析。
- `11_compare_with_GSE294965`：按 `01` → `02` → `03` / `04` 运行，生成公共数据伪 bulk 后与本项目结果比较。`01_create_pseudo_bulk` 实际写入的目录名为 `01_viz_data`，下游读取路径与此一致。

## 当前复现限制

- `07_scTM` 和 `08_gene_list` 读取 `results/06_spahdmap/` 下的预处理结果，对应上游代码未收录，运行前需要另行准备。
- `07_scTM/01_run_sctm.ipynb` 中保存 `scTM_15topics.h5ad` 的语句被注释，而后续主题比较会读取该文件；需确认对象内容并启用保存，或提供已有结果。
- `08_gene_list/01_viz_gene_score.ipynb` 将样本列表再次设置为仅 `X1912`。生成全部样本结果前需要调整该设置。
- 部分人工标注 Notebook 保留了针对特定样本的图层名和输出名；例如 `03_identify_podocyte_niche` 读取 `X1845`，但保存文件仍命名为 `21022`，复用前应核对。
- Notebook 编号表示原始分析组织方式，分支间仍需按实际输入、输出关系运行。已有绘图或执行输出不代表能在新环境直接复现。

## 主要输出

结果主要写入仓库父目录的 `results/<分析目录>/<步骤>/`：

- `.h5ad`：质控后的空间表达、参考图谱、细胞丰度、聚类及伪 bulk 对象。
- `.csv` / `.tsv`：计数矩阵、样本元数据、差异表达、细胞组成、空间簇标签和主题权重。
- `.zarr`：用于交互式组织图像与区域标注的 SpatialData 对象。
- `.png` / `.pdf`：空间分布、PCA/UMAP、热图、火山图和组间比较图；部分图仅显示在 Notebook 输出中。

## 许可证

代码采用 [MIT License](LICENSE)，Copyright © 2026 Zhijian Li。
