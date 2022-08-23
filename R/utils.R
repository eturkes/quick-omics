# Copyright 2020, 2022 Emir Turkes, UK DRI at UCL
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# This file holds common functions and methods.
# ---------------------------------------------

#' ggplot2 function providing custom aesthetics and automatic placement of categorical labels.
#' For continuous data, a colorbar is implemented.
#'
#' @param data SingleCellExperiment or Seurat object.
#' @param x,y Dimensionality reduction coordinates.
#' @param color Column metadata to color points by.
#' @param type \code{"cat"} is categorical, \code{"cont"} is continuous, \code{"NULL"} is generic.
#' @examples
#' red_dim_plot(data = sce, x = "tsne1", y = "tsne2", color = "cluster", type = "cat")
#' red_dim_plot(data = seurat, x = "umap1", y = "umap2", color = "nUMI", type = "cont")
#'
red_dim_plot <- function(data, x, y, color, type = NULL) {

  if ((class(data))[1] == "SingleCellExperiment") {
    gg_df <- data.frame(colData(data)[ , c(x, y, color)])
  } else if ((class(data))[1] == "Seurat") {
    gg_df <- data.frame(data[[x]], data[[y]], data[[color]])
  }
  rownames(gg_df) <- NULL
  gg_df[[color]] <- factor(gg_df[[color]])

  gg <- ggplot(gg_df, aes_string(x = x, y = y, color = color)) +
    geom_point(alpha = 0.35, stroke = 0.05, shape = 21, aes_string(fill = color)) +
    theme_classic() +
    theme(
      legend.position = "right", plot.title = element_text(hjust = 0.5),
      legend.title = element_blank()
    ) +
    guides(color = guide_legend(override.aes = list(alpha = 1)))

  if (is.null(type)) {
    return(gg)

  } else if (type == "cat") {
    label_df <- gg_df %>% group_by_at(color) %>% summarise_at(vars(x:y), median)
    label_df <- cbind(label_df[[1]], label_df)
    names(label_df) <- c("label", color, x, y)
    gg <- gg + geom_label_repel(data = label_df, aes(label = label), show.legend = FALSE)

  } else if (type == "cont") {
    # TODO: Refactor repeated code.
    if ((class(data))[1] == "SingleCellExperiment") {
      gg_df <- data.frame(colData(data)[ , c(x, y, color)])
    } else if ((class(data))[1] == "Seurat") {
      gg_df <- data.frame(data[[x]], data[[y]], data[[color]])
    }
    rownames(gg_df) <- NULL

    gg <- ggplot(gg_df, aes_string(x = x, y = y)) +
      geom_point(alpha = 0.35, stroke = 0.05, aes_string(color = color)) +
      theme_classic() +
      theme(
        legend.position = "right", plot.title = element_text(hjust = 0.5),
        legend.title = element_blank()
      ) +
      scale_color_viridis()
  }
  gg
}

#' Adds download buttons and horizontal scrolling to `DT::datatable`
#'
#' @param dt A data.table object.
#' @examples
#' datatable_custom(dt = data_table)
#'
datatable_custom <- function(dt) {

  datatable(
    dt, extensions = "Buttons",
    options = list(
      scrollX = TRUE, dom = "Blfrtip",
      buttons = list(
        "copy", "print",
        list(extend = "collection", buttons = c("csv", "excel", "pdf"), text = "Download")
      )
    )
  )
}

#' Convert human to mouse gene names.
#' Adapted from:
#' https://rjbioinformatics.com/2016/10/14/converting-mouse-to-human-gene-names-with-biomart-package/
#'
#' @param genes A list of human genes.
#' @examples
#' human_to_mouse_genes(genes = gene_list)
#'
human_to_mouse_genes <- function(genes) {

  human <- useMart(
    "ensembl", dataset = "hsapiens_gene_ensembl", host = "http://useast.ensembl.org/"
  )
  mouse <- useMart(
    "ensembl", dataset = "mmusculus_gene_ensembl", host = "http://useast.ensembl.org/"
  )

  new_genes <- getLDS(
    attributes = "hgnc_symbol", filters = "hgnc_symbol", values = genes , mart = human,
    attributesL = "mgi_symbol", martL = mouse, uniqueRows = TRUE
  )
  new_genes <- unique(new_genes[ , 2])
}

#' Pipeline for clustering and dimensionality reduction of post-QC scRNA-seq data.
#'
#' @param seurat Post-QC Seurat object.
#' @param assets_dir Directory to save post-processed Seurat object.
#' @param analysis_no Integer to save object to a specific directory within \code{"assets_dir"}.
#' @param sub_name Subset level for naming of cache object.
#' @param organism \code{"human"} for human, \code{"mouse"} for mouse.
#' @param vars_to_regress Vectors of nuisance variables to regress out.
#' @param parallel_override See function \code{"parallel_plan"}.
#' @examples
#' cluster_pipeline(
#'   seurat = seurat, assets_dir = assets_dir, analysis_no = 1, sub_name = "neuronal",
#'   organism = "mouse", vars_to_regress = c(batch, "mito_percent"), parallel_override = NULL)
#' )
cluster_pipeline <- function(
  seurat, assets_dir, analysis_no, sub_name, organism, vars_to_regress, parallel_override
) {

  # Utilize a basic caching system.
  rds <- file.path(assets_dir, "cache", paste0("0", analysis_no), paste0(sub_name, "_seurat.rds"))
  if (file.exists(rds)) {
    seurat <- readRDS(rds)
    return(seurat)

  } else {
    # Perform sctransform.
    # Note that this function produces many iterations of the following benign warning:
    # Warning in theta.ml(y = y, mu = fit$fitted): iteration limit reached
    DefaultAssay(seurat) <- "RNA" # Calculate off raw data, not previous "SCT" slot.
    parallel_plan(seurat, parallel_override = parallel_override)
    seurat <- SCTransform(seurat, vars.to.regress = vars_to_regress, verbose = FALSE)

    # Perform cell cycle scoring.
    if (organism == "human") {
      s_genes <- cc.genes.updated.2019$s.genes
      g2m_genes <- cc.genes.updated.2019$g2m.genes
    } else if (organism == "mouse") {
      s_genes <- human_to_mouse_genes(cc.genes.updated.2019$s.genes)
      g2m_genes <- human_to_mouse_genes(cc.genes.updated.2019$g2m.genes)
    }
    seurat <- CellCycleScoring(seurat, s.features = s_genes, g2m.features = g2m_genes)
    seurat$cc_diff <- seurat$S.Score - seurat$G2M.Score # Combined proliferating cell signal.

    # Perform PCA.
    # We also add the output as column metadata for use with `red_dim_plot`.
    seurat <- RunPCA(seurat, verbose = FALSE)
    add_df <- data.frame(Embeddings(seurat, reduction = "pca")[ , 1:2])
    names(add_df) <- paste0("pca", seq(ncol(add_df)))
    seurat$pca1 <- add_df$pca1
    seurat$pca2 <- add_df$pca2

    # Perform UMAP reduction.
    # We also add the output as column metadata for use with `red_dim_plot`.
    seurat <- RunUMAP(seurat, dims = 1:30, min.dist = 0.75, verbose = FALSE)
    add_df <- data.frame(Embeddings(seurat, reduction = "umap"))
    names(add_df) <- paste0("umap", seq(ncol(add_df)))
    seurat$umap1 <- add_df$umap1
    seurat$umap2 <- add_df$umap2

    # Perform Louvain clustering.
    resolution <- (dim(seurat)[2] / 3000) * 0.8 # Default is optimal for 3K cells so we scale it.
    seurat <- FindNeighbors(seurat, dims = 1:30, verbose = FALSE)
    seurat <- FindClusters(seurat, resolution = resolution, verbose = FALSE)

    saveRDS(seurat, rds)
  }
  seurat
}

#' Set the `plan` for `future` based on free memory and object size with the option to override.
#'
#' @param object Object to check if \code{"future.globals.maxSize"} large enough to parallelize.
#' @param parallel_override \code{"NULL"} to calculate plan decision, \code{0} for sequential, a
#' non-zero integer for multiprocess and to set `future.globals.maxSize`.
#' @examples
#' parallel_plan(object = seurat, parallel_override = 5368709120)
#'
parallel_plan <- function(object, parallel_override = NULL) {

  if (is.null(parallel_override)) {
    # Get free memory.
    gc()
    mem <- as.numeric(unlist(strsplit(system("free -b", intern = TRUE)[2], " "))[7])

    # Distribute free memory (minus 10 GiB) across available cores.
    mem <- mem - 10 * 1024 ^ 3
    mem <- mem / detectCores()

    # Enable parallelization only if `object` can fit in `future.globals.maxSize` (plus 1 Gib).
    if (mem > object.size(object) + 1 * 1024 ^ 3) {
      plan("multiprocess")
      options(future.globals.maxSize = mem)
    } else {
      plan("sequential")
    }

  } else if (parallel_override == 0) {
    plan("sequential")

  } else {
    plan("multiprocess")
    options(future.globals.maxSize = parallel_override)
  }
}
