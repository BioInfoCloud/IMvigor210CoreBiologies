# IMvigor210CoreBiologies 修复指南（针对过时的 DESeq 依赖）

> 适用版本：IMvigor210CoreBiologies 1.0.0（2019-02-11 发布，原始依赖已大面积失效）
> 状态：**已实际修复并验证**（`R CMD INSTALL` 成功，`data(cds)` 返回可用的 `DESeqDataSet`，`counts/pData/fData` 与 DESeq2 工作流均可用）。

## 1. 问题诊断

这个包之所以「装不上 / 装上了也用不了」，根因有两个，且**代码本身其实并不依赖老版 DESeq**：

### 1.1 `DESCRIPTION` 里声明的 `DESeq` 已被 Bioconductor 移除
- 老的 `DESeq` 包（与 `DESeq2` 不同）在 **Bioconductor 3.13** 被正式 removed，新版本 R 已无法从仓库安装。
- 它出现在 `Imports:` 字段里，导致 `install.packages(..., repos=NULL)` / `R CMD INSTALL` 一开始就因「缺少 DESeq 依赖」而失败。
- **但包里的 R 函数并没有真正调用老 `DESeq`**：`R/customFunctions.R` 的 `mySimpleVoom()` 虽然有个形参叫 `DESeq`，实际只调用 `DESeq2::estimateSizeFactorsForMatrix()`、`edgeR::*` 和 `limma::voom()`。Grep 全局确认：没有任何 `DESeq::` 命名空间调用。

### 1.2 真正的死依赖是数据对象 `cds`
- 包随附的 `data/cds.RData` 把表达数据存成了老 `DESeq` 的 **`CountDataSet`** S4 类。
- 关键发现：`CountDataSet` 的插槽（`assayData` / `phenoData@data` / `featureData@data`）**即使不加载 DESeq 命名空间也能通过 S4 插槽访问**（会有一条 benign warning：`namespace 'DESeq' is not available and has been replaced by .GlobalEnv`，可忽略）。因此我们可以**完全不安装那个已死的 DESeq 包**就把数据抽取并重建为现代对象。

### 1.3 另一个已归档的包：`lsmeans`
- `lsmeans` 已被 CRAN 归档，被 `emmeans` 取代。它只用于 `inst/analysis/` 下的小鼠影像脚本，并不被包的导出函数使用，但同样列在 `Imports` 里挡住了安装。

## 2. 本次已完成的修改

| 文件 | 修改 |
|------|------|
| `DESCRIPTION` | 从 `Imports` 删除 `DESeq`；把只用于分析脚本的 `lsmeans / circlize / ComplexHeatmap / corrplot / DT / plyr / reshape2 / spatstat` 移至 `Suggests`；为支持 re-export 的访问器，新增 `SummarizedExperiment, S4Vectors` 到 `Imports`。保留 `biomaRt, DESeq2, dplyr, edgeR, ggplot2, graphics, limma, methods, stats, survival`。 |
| `inst/analysis/imageAnalysis/mouseTumorImaging-{1430,1436,666}-CD3-dataManipulation.r` | `library(lsmeans)` → `library(emmeans)` |
| `inst/analysis/mouseTumorImaging-CD3-pooledAnalysis.Rmd` | `library(lsmeans)` → `library(emmeans)`；`lsmeans(fit.sum, ~ Treatment)` → `emmeans(...)` |
| `vignettes/index.Rmd` | 安装说明移除 `DESeq`；废弃的 `biocLite()` 改为 `BiocManager::install()`；`lsmeans` → `emmeans`；`CountDataSet` 文案改为 `DESeqDataSet` |
| `R/data.R`（及生成的 `man/cds.Rd`） | 文档文案 `countDataSet` → `DESeqDataSet` |
| `R/reexports.R`（新增） | 把 `counts`（DESeq2）、`colData/rowData/assay/assays`（SummarizedExperiment）、`metadata`（S4Vectors）re-export，使用户在 `library()` 后无需再手动 `library(DESeq2)` 即可直接 `counts(cds)` 等 |
| `R/compatibility.R`（新增） | `.onLoad` 时自动为 `DESeqDataSet` 注册 `pData`/`fData` 方法（分别映射到 `colData`/`rowData`），历史代码 `pData(cds)`/`fData(cds)` 无需任何额外操作即可用 |
| `inst/migrate_cds_to_DESeq2.R`（重写） | **不再需要安装老 DESeq**：直接通过 S4 插槽抽取 `counts/pData/fData`，重建为 `DESeqDataSet` 并写回 `data/cds.RData`（注意 `save(cds = dds, ...)` 会错误地以 `dds` 为名存储，脚本已改用 `cds <- dds; save(cds, ...)`）。原始 `CountDataSet` 已备份在包外 `F:/DatabaseData/IMvigor210/cds.RData.orig.CountDataSet_backup`，避免被 `R CMD build` 一并打包。 |
| `inst/compatibility_shim.R`（增强） | 独立可用的兼容层（方法签名已修正为与 Biobase 泛型 `(object)` 一致）。现已被 `.onLoad` 自动覆盖，作为离线/手动备用保留。 |

> 完成以上修改后，`R CMD INSTALL` 已成功把包装上；`data(cds)` 现返回现代 `DESeqDataSet`，所有访问器与 DESeq2 工作流均可用。

## 3. 迁移 `cds` 数据对象（关键，且已实际执行）

脚本 `inst/migrate_cds_to_DESeq2.R` 已重写并**实际跑通**（在本地 R 4.6.1 + 已安装的 Bioconductor 环境下，把 `data/cds.RData` 由 `CountDataSet` 迁移为 `DESeqDataSet`）。

### 3.1 前置依赖（全部是现代、可正常安装的包）
```r
BiocManager::install(c("DESeq2", "SummarizedExperiment", "S4Vectors", "Biobase"))
```
> **不再需要安装任何已废弃的包**（原来的老 `DESeq` 完全不需要）。

### 3.2 运行迁移（在包根目录执行）
```r
Rscript inst/migrate_cds_to_DESeq2.R
```
脚本会：
1. 自动把原 `data/cds.RData` 备份为 `data/cds.RData.orig`；
2. 通过 S4 插槽读取旧 `CountDataSet`，抽取 `counts / pData / fData`（无需 DESeq）；
3. 重建为现代 **`DESeqDataSet`** 并以正确的名字 `cds` 覆盖写回 `data/cds.RData`；
4. 新对象的访问方式：`counts(cds)`（不变）、`colData(cds)`（原 `pData`）、`rowData(cds)`（原 `fData`）。

### 3.3 兼容旧代码（现已自动生效）
历史分析脚本里若仍调用 `pData(cds)` / `fData(cds)`，在 `library(IMvigor210CoreBiologies)` 后**已自动可用**（由 `.onLoad` 注册）。如需在包作用域之外手动启用，也可：
```r
source(system.file("compatibility_shim.R", package = "IMvigor210CoreBiologies"))
```

## 4. 网上有没有人修改过这个包？

结论：**有，但都不是「干净移除 DESeq」的现代改写**，多数只是把老 `DESeq` 硬塞进去。

1. **GitHub fork：`SiYangming/IMvigor210CoreBiologies`**（2021-02）
   - 做法：把 `DESeq_1.39.0` 的预编译包挂到 GitHub Release，安装时手动本地装这个 `DESeq`，再装包。
   - 评价：**并没有删除 DESeq 依赖**，而是把死包托管到了 GitHub，且绑死了特定 R/Bioconductor 版本，换环境照样容易失败。不推荐作为长期方案。

2. **社区主流「笨办法」**（CSDN / 简书 / 腾讯云大量博客）
   - 手动从 Bioconductor 归档下载 `DESeq_1.38.0 / 1.39.0` 的 `.tar.gz`，`install.packages(..., repos=NULL)` 本地安装，再装本包。
   - 评价：能用，但脆弱、版本强耦合，且本质还是在用已废弃的包。

3. **最省事的数据获取法：不装包，直接读 `.RData`**
   - 很多教程指出：你真正要的只是数据。包源码 `data/` 目录里就是现成的 `cds.RData`，可直接 `load()`，无需安装包：
     ```r
     e <- new.env(); load("data/cds.RData", e); cds <- e$cds
     # 借助插槽直接取数（无需 DESeq）：
     counts_mat <- cds@assayData$counts
     pd  <- cds@phenoData@data
     fd  <- cds@featureData@data
     ```
   - 本修复方案的迁移脚本正是基于这一思路，但进一步把它固化成了标准 `DESeqDataSet`。

4. **官方同源的现代替代品：`easierData` 包（Bioconductor）**
   - `easierData::get_Mariathasan2018_PDL1_treatment()` 返回的是标准的 **`SummarizedExperiment`**，数据同源（同一篇 Nature 2018 论文），但样本注释略有删减（并非 25 列全量）。
   - 适合「只想要数据、想用现代对象」的场景：
     ```r
     BiocManager::install("easierData")
     library(easierData); library(SummarizedExperiment)
     se <- get_Mariathasan2018_PDL1_treatment()
     assay(se)                       # 表达矩阵
     as.data.frame(colData(se))      # 样本注释
     ```

## 5. 推荐落地方案（按需求选）

| 你的目标 | 推荐做法 |
|----------|----------|
| 想让原包 `library(IMvigor210CoreBiologies)` + `data(cds)` 正常工作 | 采用本文第 2、3 节（已完成并验证）：改 `DESCRIPTION` + 跑一次 `migrate_cds_to_DESeq2.R` |
| 只想要表达矩阵 + 临床注释做自己的分析 | 直接 `load("data/cds.RData")` 用插槽取数，或改用 `easierData` |
| 想跑 `inst/analysis` 里的小鼠影像脚本 | 已把 `lsmeans` 换成 `emmeans`，确保 `emmeans` 已装即可 |
| 想要最省心、零改造 | 用 `easierData`（现代 `SummarizedExperiment`），放弃原包 |

## 6. 验证（本机已实测通过）

```r
library(IMvigor210CoreBiologies)
data(cds)
class(cds)                          # "DESeqDataSet"
dim(counts(cds))                    # 31286 features x 348 samples
dim(colData(cds))                   # 348 x 25   （原 pData）
dim(rowData(cds))                   # 31286 x 6  （原 fData）
pData(cds)                          # 自动可用（.onLoad 注册），348 x 25 data.frame
fData(cds)                          # 自动可用，31286 x 6 data.frame

# 现代 DESeq2 工作流可直接跑：
suppressMessages(library(DESeq2))
dds <- estimateSizeFactors(cds)     # 348 个 size factors 全部有限
```

实测结果：`R CMD INSTALL` 成功；`data(cds)` 为 `DESeqDataSet`；`counts`/`colData`/`rowData` 经 re-export 直接可用；`pData`/`fData` 经 `.onLoad` 自动可用；`estimateSizeFactors()` 正常运行。

## 7. 注意事项
- `NAMESPACE` 已通过 `roxygen2::roxygenise()` 重新生成（新增了 re-export 与 `.onLoad` 相关条目）。若后续再改 R/ 文档，记得重新 `roxygenise()`。
- 文档文案（`R/data.R`、`man/cds.Rd`、`vignettes/index.Rmd`）已从 `CountDataSet` 更正为 `DESeqDataSet`。
- 原始 `CountDataSet` 备份位于包外 `F:/DatabaseData/IMvigor210/cds.RData.orig.CountDataSet_backup`（51MB），建议保留以便回滚；它已移出包目录，不会被 `R CMD build` 打包。
- 若你只想临时读数据、不想动包：`load("data/cds.RData")` + 插槽取数即可，完全不需要 DESeq。
