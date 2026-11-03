# MatriCom — input readers (RDS / QS / H5AD).

# --- uploads ---------------------------------------------------------------

read_qs_any <- function(path) {
  # qs2 format first; legacy 'qs' files only if the archived qs package is
  # present (it is no longer on CRAN, so it is not a declared dependency).
  out <- tryCatch(qs2::qs_read(path), error = function(e) e)
  if (!inherits(out, "error")) return(out)
  legacy <- "qs"   # held in a variable so deployment tools don't require it
  if (requireNamespace(legacy, quietly = TRUE)) return(getExportedValue(legacy, "qread")(path))
  stop("This .qs file was written with the retired 'qs' package. Re-save it with ",
       "qs2::qs_save() or saveRDS() and upload again.")
}

h5_attr <- function(node, name, default = NULL) {
  if (name %in% hdf5r::h5attr_names(node)) hdf5r::h5attr(node, name) else default
}

# read one anndata element: plain dataset, categorical (codes/categories),
# nullable array (values/mask, anndata >= 0.11) or legacy __categories column
h5_read_elem <- function(grp, nm) {
  node <- grp[[nm]]
  if (inherits(node, "H5Group")) {
    kids <- names(node)
    if (all(c("codes", "categories") %in% kids)) {
      codes <- as.integer(node[["codes"]]$read())
      cats <- as.character(h5_read_elem(node, "categories"))
      v <- rep(NA_character_, length(codes)); ok <- !is.na(codes) & codes >= 0
      v[ok] <- cats[codes[ok] + 1L]
      return(factor(v, levels = cats))
    }
    if (all(c("values", "mask") %in% kids)) {
      v <- node[["values"]]$read()
      mask <- as.logical(node[["mask"]]$read())
      v[mask] <- NA
      return(v)
    }
    return(NULL)
  }
  v <- node$read()
  if (length(dim(v)) > 1) return(NULL)
  if ("__categories" %in% names(grp) && nm %in% names(grp[["__categories"]])) { # anndata 0.7
    cats <- as.character(grp[["__categories"]][[nm]]$read())
    v2 <- rep(NA_character_, length(v)); ok <- !is.na(v) & v >= 0
    v2[ok] <- cats[v[ok] + 1L]
    return(factor(v2, levels = cats))
  }
  v
}

h5_index <- function(grp) {
  idx <- h5_attr(grp, "_index", NULL)
  cand <- unique(c(idx, "_index", "index", "obs_names", "var_names"))
  for (nm in cand) if (!is.null(nm) && nm %in% names(grp)) return(as.character(h5_read_elem(grp, nm)))
  NULL
}

h5_frame <- function(node) {
  if (inherits(node, "H5D")) {                                  # anndata <= 0.6 compound
    df <- as.data.frame(node$read())
    ix <- intersect(c("index", "_index"), names(df))
    if (length(ix)) { rownames(df) <- as.character(df[[ix[1]]]); df[[ix[1]]] <- NULL }
    return(df)
  }
  idx <- h5_index(node)
  ord <- h5_attr(node, "column-order", NULL)
  cols <- setdiff(names(node), c("__categories", h5_attr(node, "_index", "_index"), "_index"))
  if (!is.null(ord) && length(ord)) cols <- c(intersect(as.character(ord), cols), setdiff(cols, ord))
  df <- data.frame(row.names = if (is.null(idx)) NULL else make.unique(idx))
  n <- length(idx)
  for (nm in cols) {
    v <- tryCatch(h5_read_elem(node, nm), error = function(e) NULL)
    if (!is.null(v) && (is.null(idx) || length(v) == n)) df[[nm]] <- v
  }
  df
}

h5_matrix <- function(node, n_obs, n_var) {
  # returns var x obs (genes x cells)
  if (inherits(node, "H5Group")) {
    enc <- h5_attr(node, "encoding-type", h5_attr(node, "h5sparse_format", "csr_matrix"))
    shape <- h5_attr(node, "shape", h5_attr(node, "h5sparse_shape", c(n_obs, n_var)))
    x <- as.numeric(node[["data"]]$read())
    i <- as.integer(node[["indices"]]$read())
    p <- as.integer(node[["indptr"]]$read())
    if (grepl("csc", enc)) {
      m <- Matrix::sparseMatrix(i = i + 1L, p = p, x = x, dims = c(shape[1], shape[2]))
      return(Matrix::t(m))
    }
    return(Matrix::sparseMatrix(i = i + 1L, p = p, x = x, dims = c(shape[2], shape[1])))
  }
  m <- node$read()        # hdf5r returns HDF5 (obs, var) as an R var x obs array
  if (!identical(dim(m), c(as.integer(n_var), as.integer(n_obs)))) m <- t(m)
  as(Matrix::Matrix(m, sparse = TRUE), "CsparseMatrix")
}

read_h5ad <- function(path) {
  h5 <- hdf5r::H5File$new(path, mode = "r")
  on.exit(h5$close_all())
  obs <- h5_frame(h5[["obs"]])
  var <- h5_frame(h5[["var"]])
  n_obs <- nrow(obs); n_var <- nrow(var)
  # raw counts preferred: layers/counts, then raw/X, then X
  if ("layers" %in% names(h5) && "counts" %in% names(h5[["layers"]])) {
    X <- h5_matrix(h5[["layers"]][["counts"]], n_obs, n_var); genes <- rownames(var)
  } else if ("raw" %in% names(h5) && "X" %in% names(h5[["raw"]])) {
    rv <- h5_frame(h5[["raw"]][["var"]])
    X <- h5_matrix(h5[["raw"]][["X"]], n_obs, nrow(rv)); genes <- rownames(rv)
  } else {
    X <- h5_matrix(h5[["X"]], n_obs, n_var); genes <- rownames(var)
  }
  dimnames(X) <- list(genes, rownames(obs))
  build_from_counts(X, obs)
}

# Any supported upload -> list(data = genes x cells normalised matrix,
#                              meta = cell metadata, assays = assay names)
load_input <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "h5ad") return(read_h5ad(path))
  obj <- if (ext == "qs") read_qs_any(path) else if (ext == "rds") readRDS(path) else stop("Unsupported file type.")

  if (inherits(obj, "Seurat")) {
    assays <- names(obj@assays)
    if (!"RNA" %in% assays) return(list(data = NULL, meta = obj@meta.data, assays = assays))
    major <- strsplit(as.character(obj@version), split = "\\.")[[1]][1]
    if (major == "3") {
      return(list(data = obj@assays$RNA@data, meta = obj@meta.data, assays = assays))
    }
    rna <- obj@assays[["RNA"]]
    if (inherits(rna, "Assay5") && length(SeuratObject::Layers(rna, search = "counts")) > 1)
      rna <- SeuratObject::JoinLayers(rna)
    cnt <- SeuratObject::LayerData(rna, layer = "counts")
    return(build_from_counts(cnt, obj@meta.data))
  }
  if (inherits(obj, "SingleCellExperiment")) {
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) stop("SingleCellExperiment not installed")
    cnt <- SummarizedExperiment::assay(obj, "counts")
    meta <- as.data.frame(SummarizedExperiment::colData(obj))
    return(build_from_counts(cnt, meta))
  }
  if (is.list(obj) && all(c("X", "obs", "var") %in% names(obj))) {
    cnt <- obj$X
    f <- obj$var
    if (!is.null(f) && "gene_ids" %in% colnames(f)) rownames(cnt) <- f$gene_ids
    else if (!is.null(f) && "index" %in% colnames(f)) rownames(cnt) <- f$index
    return(build_from_counts(cnt, obj$obs))
  }
  stop("Unknown object structure.")
}
