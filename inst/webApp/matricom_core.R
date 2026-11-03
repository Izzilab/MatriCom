# MatriCom core — vectorised re-implementation of matricom.slim(), filtres(),
# matrienrich() and the influencer computation of the original app.
#
# Output columns and interaction calls follow the original MatriCom; the
# corrections made in the 2026 rebuild are listed in CHANGES.md.

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})

NO_INT_MSG <- "No interactions found with these parameters. Please check your selection!"

fail_df <- function(msg = NO_INT_MSG) {
  df <- data.frame(error.message = msg)
  attributes(df)$outcome <- "failure"
  df
}

is_failure <- function(x) {
  o <- attributes(x)$outcome
  !is.null(o) && identical(o, "failure")
}

# ---------------------------------------------------------------------------
# Input normalisation helpers
# ---------------------------------------------------------------------------

as_dgC <- function(m) {
  if (inherits(m, "dgCMatrix")) return(m)
  if (inherits(m, "DelayedArray") || inherits(m, "DelayedMatrix")) {
    m <- as(m, "dgCMatrix")
    return(m)
  }
  as(as(as(m, "CsparseMatrix"), "generalMatrix"), "dMatrix")
}

# Seurat::LogNormalize (scale.factor = 1e4), identical operation order:
# log1p(x / colSum * 1e4)
log_normalize <- function(counts) {
  counts <- as_dgC(counts)
  cs <- colSums(counts)
  nn <- diff(counts@p)
  counts@x <- log1p(counts@x / rep.int(cs, nn) * 10000)
  counts
}

# Equivalent of the original upload path for non-v3 objects:
#   CreateAssayObject(counts) -> CreateSeuratObject(meta.data) -> NormalizeData
# SeuratObject (not Seurat) is used so feature/cell-name sanitisation and the
# automatic metadata columns (orig.ident, nCount_RNA, nFeature_RNA) match.
build_from_counts <- function(counts, meta = NULL) {
  old <- options(Seurat.object.assay.version = "v3"); on.exit(options(old))
  counts <- as_dgC(counts)
  assay <- suppressWarnings(SeuratObject::CreateAssayObject(counts = counts))
  so <- suppressWarnings(SeuratObject::CreateSeuratObject(counts = assay, meta.data = meta))
  cnt <- SeuratObject::LayerData(so[["RNA"]], layer = "counts")
  list(data = log_normalize(cnt), meta = so@meta.data, assays = "RNA")
}

# Ensembl -> symbol conversion (EnsDb.Hsapiens.v86 lookup, precomputed into
# ens2sym.rds so the 60 MB annotation package is not needed at run time).
# Mirrors the original: unmapped IDs get one shared random placeholder.
convert_ensembl <- function(ids, ens2sym) {
  sym <- ens2sym$SYMBOL[match(ids, ens2sym$GENEID)]
  sym[is.na(sym)] <- paste0("not.annotated.gene_", sample(c(1:1000000), 1))
  sym
}

# ---------------------------------------------------------------------------
# matricom.slim()
#   mat     : genes x cells normalised matrix (the RNA "data" slot)
#   groups  : per-cell identity values (metadata column), same order as columns
# ---------------------------------------------------------------------------

matricom_core <- function(mat, groups, min.pct, expr.filter, lr, mlist,
                          homo.policy = FALSE, conv = FALSE, ens2sym = NULL) {

  mat <- as_dgC(mat)
  if (is.null(colnames(mat)) || length(colnames(mat)) < 1)
    colnames(mat) <- paste0("sample.", seq_len(ncol(mat)))
  if (is.null(rownames(mat)) || length(rownames(mat)) < 1)
    rownames(mat) <- paste0("gene.", seq_len(nrow(mat)))

  rownames(mat) <- toupper(rownames(mat))          # "enforcing humanity"
  if (isTRUE(conv)) rownames(mat) <- convert_ensembl(rownames(mat), ens2sym)

  lr_genes <- unique(c(lr$V1, lr$V2))

  # subset(features = LR genes, idents = non-NA groups); first match on dups
  feats <- intersect(lr_genes, rownames(mat))
  if (length(feats) == 0) stop("Cannot find features provided")
  keep <- !is.na(groups)
  mat <- mat[match(feats, rownames(mat)), keep, drop = FALSE]
  g <- make.names(gsub("-", "_", as.character(groups[keep])))

  # drop populations with a single cell
  tab <- table(g)
  big <- names(tab)[tab > 1]
  if (length(big) < length(tab)) {
    k <- g %in% big
    mat <- mat[, k, drop = FALSE]
    g <- g[k]
  }
  if (ncol(mat) == 0) stop("no cells")

  pops <- unique(g)
  gi <- match(g, pops)
  ncell <- tabulate(gi, nbins = length(pops))

  # one pass over the sparse matrix: per gene x population
  #   npos  = # cells with value > 0      (get_expressed_genes / Percent_Expressing)
  #   nnz   = # cells with value != 0     (mean of non-zero values)
  #   ssum  = sum of values
  T <- as(mat, "TsparseMatrix")
  ri <- T@i + 1L
  pj <- gi[T@j + 1L]
  key <- (pj - 1L) * nrow(mat) + ri
  dims <- c(nrow(mat), length(pops))
  npos <- matrix(tabulate(key[T@x > 0], nbins = prod(dims)), dims[1], dims[2])
  nz <- T@x != 0
  nnz  <- matrix(tabulate(key[nz], nbins = prod(dims)), dims[1], dims[2])
  ssum <- matrix(0, dims[1], dims[2])
  if (any(nz)) {
    agg <- rowsum(T@x[nz], key[nz], reorder = TRUE)
    ssum[as.integer(rownames(agg))] <- agg[, 1]
  }
  rm(T)
  dimnames(npos) <- dimnames(nnz) <- dimnames(ssum) <- list(feats, pops)

  frac <- sweep(npos, 2, ncell, "/")            # sum(x > 0) / n
  pct <- min.pct / 100

  ligavals <- list()
  for (p in seq_along(pops)) {
    expressed <- frac[, p] >= pct
    if (sum(expressed) < 2) next                 # original: nrow(k) < 2 -> next
    mn <- ssum[expressed, p] / nnz[expressed, p] # mean of non-zero values
    mn <- mn[mn >= expr.filter]
    if (length(mn) < 1) next
    ligavals[[pops[p]]] <- mn
  }
  ligands <- lapply(ligavals, function(v) toupper(names(v)))
  if (length(ligands) == 0) return(fail_df())

  # fraction of cells with value > 0 (same quantity used to call a gene expressed)
  pexp <- frac

  is_m1 <- lr$V1 %in% mlist$gene
  is_m2 <- lr$V2 %in% mlist$gene
  lr_ok <- is_m1 | is_m2
  lr_type <- ifelse(is_m1 & is_m2, "extracellular", "cell-matrix")

  # membership of every LR gene in every population's ligand set
  lp <- names(ligands)
  M1 <- vapply(ligands, function(s) lr$V1 %in% s, logical(nrow(lr)))
  M2 <- vapply(ligands, function(s) lr$V2 %in% s, logical(nrow(lr)))
  if (!is.matrix(M1)) { M1 <- matrix(M1, ncol = length(lp)); M2 <- matrix(M2, ncol = length(lp)) }
  M1[!lr_ok, ] <- FALSE
  M2[!lr_ok, ] <- FALSE

  # ---- homocellular (sender == receiver), population order then LR order
  homo_idx <- which(M1 & M2, arr.ind = TRUE)
  homo_idx <- homo_idx[order(homo_idx[, 2], homo_idx[, 1]), , drop = FALSE]
  omo <- data.table(V1 = lr$V1[homo_idx[, 1]], V2 = lr$V2[homo_idx[, 1]],
                    source = lr$source[homo_idx[, 1]], relscore = lr$relscore[homo_idx[, 1]],
                    final = lr_type[homo_idx[, 1]],
                    sender = lp[homo_idx[, 2]], receiver = lp[homo_idx[, 2]])
  if (nrow(omo) > 0) {
    # original: ints_omo[!duplicated(t(apply(ints_omo, 1, sort))), ]
    # (duplicate = same multiset of the 7 cell values, numbers via format())
    cm <- cbind(omo$V1, omo$V2, omo$source, format(omo$relscore), omo$final,
                omo$sender, omo$receiver)
    long <- data.table(r = rep(seq_len(nrow(cm)), 7L), v = as.vector(cm))
    setorder(long, r, v)
    keys <- long[, .(k = paste(v, collapse = "\r")), by = r]$k
    omo <- omo[!duplicated(keys)]
  }

  # ---- heterocellular: for sender i, receiver w (i != w), LR rows with V1 in
  # L_i and V2 in L_w. (The original also built the mirrored et_2 block; every
  # et_2 row that survives its later filters is an exact duplicate of an et_1
  # row of the same block and is removed by the final distinct().)
  npop <- length(lp)
  het <- NULL
  if (npop > 1) {
    a <- which(M1, arr.ind = TRUE)  # (lr row, sender)
    b <- which(M2, arr.ind = TRUE)  # (lr row, receiver)
    A <- data.table(r = a[, 1], i = a[, 2])
    B <- data.table(r = b[, 1], w = b[, 2])
    het <- A[B, on = "r", allow.cartesian = TRUE, nomatch = NULL]
    het <- het[i != w]
    setorder(het, i, w, r)
  }
  et <- if (!is.null(het) && nrow(het) > 0) {
    data.table(V1 = lr$V1[het$r], V2 = lr$V2[het$r], source = lr$source[het$r],
               relscore = lr$relscore[het$r], final = lr_type[het$r],
               sender = lp[het$i], receiver = lp[het$w])
  } else NULL

  fin <- rbindlist(list(omo, et))
  if (nrow(fin) < 1) return(fail_df())
  setnames(fin, c("Gene1", "Gene2", "source", "relscore", "Type of interaction",
                  "Population1", "Population2"))
  setcolorder(fin, c(6, 1, 5, 2, 7, 3, 4))

  if (shiny::isTruthy(homo.policy)) {
    fin <- fin[!(Population1 == Population2 & Gene1 == Gene2)]
  }
  fin <- na.omit(fin)
  if (nrow(fin) < 1) return(fail_df())

  # per-row lookups (the original apply(fin, 1, ...) block)
  relscore_chr <- format(fin$relscore)            # as.matrix() formatting
  gr1 <- match(fin$Gene1, feats); pc1 <- match(fin$Population1, pops)
  gr2 <- match(fin$Gene2, feats); pc2 <- match(fin$Population2, pops)
  a1 <- pexp[cbind(gr1, pc1)]; a1[is.na(a1)] <- 0
  a2 <- pexp[cbind(gr2, pc2)]; a2[is.na(a2)] <- 0
  look <- function(pop, gene) {
    out <- numeric(length(pop))
    for (p in unique(pop)) {
      ix <- which(pop == p)
      v <- ligavals[[p]]
      if (is.null(v)) next
      m <- v[match(gene[ix], names(v))]
      m[is.na(m)] <- 0
      out[ix] <- m
    }
    out
  }
  m1 <- look(fin$Population1, fin$Gene1)
  m2 <- look(fin$Population2, fin$Gene2)

  res <- data.frame(Population1 = fin$Population1,
                    Gene1 = fin$Gene1,
                    mean.expr.Gene1 = m1,
                    perc.expr.Population1 = a1,
                    `Type of interaction` = fin$`Type of interaction`,
                    Gene2 = fin$Gene2,
                    mean.expr.Gene2 = m2,
                    perc.expr.Population2 = a2,
                    Population2 = fin$Population2,
                    Reliability.score = as.numeric(relscore_chr),
                    check.names = FALSE, stringsAsFactors = FALSE)

  res <- res[order(-res$Reliability.score), ]
  res <- res[res$perc.expr.Population1 >= pct & res$perc.expr.Population2 >= pct, ]
  res <- res[res$Reliability.score > 0, ]
  res <- res[res$mean.expr.Gene1 >= expr.filter & res$mean.expr.Gene2 >= expr.filter, ]
  res <- res[!duplicated(res), ]
  rownames(res) <- NULL
  res
}

# ---------------------------------------------------------------------------
# filtres() — post-run filters
# ---------------------------------------------------------------------------

filtres_core <- function(data, expr, exc, type, hom, cc, dups, excl, CCgenes2) {
  m <- data

  m <- m[m$Reliability.score %in% expr, ]
  if (nrow(m) < 1) return(fail_df())

  if (shiny::isTruthy(hom)) {
    m <- m[m$Gene1 != m$Gene2, ]
    if (nrow(m) < 1) return(fail_df())
  }

  if (!isFALSE(exc)) {
    # Exclusion list: a multimer (multimers.txt) must be assembled inside one
    # cell, so a heterocellular pair whose two genes belong to the same complex
    # is removed. Everything else is kept.
    ex <- lapply(excl[[1]], function(x) gsub(" ", "", unlist(strsplit(x, split = ","))))
    same_complex <- logical(nrow(m))
    for (z in ex) same_complex <- same_complex | (m$Gene1 %in% z & m$Gene2 %in% z)
    m <- m[!(m$Population1 != m$Population2 & same_complex), ]
    if (nrow(m) < 1) return(fail_df())
  }

  if (length(type) == 0) return(fail_df())
  if (length(type) == 1) {
    if (type %in% "Homocellular") m <- m[m$Population1 == m$Population2, ]
    else m <- m[m$Population1 != m$Population2, ]
    if (nrow(m) < 1) return(fail_df())
  }

  m$g1g2 <- paste0(m$Gene1, "_", m$Gene2)
  if (!isFALSE(dups)) {
    # Maximize model: one row per communication. A communication is the
    # unordered pair of endpoints (population, gene), so Population A/Gene A -
    # Population B/Gene B and its reciprocal B/B - A/A are the same event.
    # Duplicates (different sources/scores, or reciprocal entries) collapse to
    # the highest reliability score; row order is otherwise preserved.
    e1 <- paste(m$Population1, m$Gene1, sep = "\r")
    e2 <- paste(m$Population2, m$Gene2, sep = "\r")
    key <- ifelse(e1 <= e2, paste(e1, e2, sep = "\n"), paste(e2, e1, sep = "\n"))
    o <- order(match(key, unique(key)), -as.numeric(m$Reliability.score))
    keep <- o[!duplicated(key[o])]
    m <- m[sort(keep), ]
  }
  rownames(m) <- NULL

  lab <- paste0(toupper(substr(CCgenes2$label, 1, 1)), tolower(substring(CCgenes2$label, 2)))
  lab[lab == "Extracellular"] <- "Extracellular (Non-matrisome)"
  sg <- unique(CCgenes2[[1]][lab %in% cc])
  fin <- m[m$Gene1 %in% sg & m$Gene2 %in% sg, ]
  if (nrow(fin) < 1) return(fail_df())
  attributes(fin)$outcome <- "success"
  fin
}

# ---------------------------------------------------------------------------
# matrienrich() — hypergeometric signature enrichment per population pair
# ---------------------------------------------------------------------------

matrienrich_core <- function(data, fams, signs) {
  l <- rbind(data.frame(standard_name = fams[[2]], gene = fams[[1]]),
             data.frame(standard_name = fams[[3]], gene = fams[[1]]),
             data.frame(standard_name = signs[[1]], gene = signs[[2]]))
  sigs <- unique(l$standard_name)
  ct <- paste0(data$Population1, "-", data$Population2)
  pairs <- unique(ct)
  if (length(pairs) == 0) return(data.frame(populations = character(), signature = character(),
                                            overlap = integer(), p.value = numeric()))
  genes <- unique(c(data$Gene1, data$Gene2, l$gene))
  # population-pair x gene incidence
  pg <- unique(data.table(p = c(ct, ct), g = c(data$Gene1, data$Gene2)))
  P <- sparseMatrix(i = match(pg$p, pairs), j = match(pg$g, genes), x = 1,
                    dims = c(length(pairs), length(genes)))
  sg <- unique(data.table(s = l$standard_name, g = l$gene))
  S <- sparseMatrix(i = match(sg$g, genes), j = match(sg$s, sigs), x = 1,
                    dims = c(length(genes), length(sigs)))
  ov <- as.matrix(P %*% S)
  nk <- colSums(S)
  ng <- rowSums(P)
  K <- matrix(nk, nrow(ov), ncol(ov), byrow = TRUE)
  N <- matrix(ng, nrow(ov), ncol(ov))
  p <- phyper(ov - 1, K, 4048 - K, N, lower.tail = FALSE)
  data.frame(populations = rep(pairs, each = length(sigs)),
             signature = rep(sigs, times = length(pairs)),
             overlap = as.integer(as.vector(t(ov))),
             p.value = as.vector(t(p)),
             stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# Network influencers (plot, on-screen table and XLSX export share this)
# ---------------------------------------------------------------------------

influencers_core <- function(out) {
  n <- unique(data.frame(V1 = out$Gene1, V2 = out$Gene2, stringsAsFactors = FALSE))
  if (nrow(n) < 1) return(NULL)
  n <- igraph::simplify(igraph::graph_from_data_frame(n))
  va <- igraph::degree(n)
  va <- va[order(-va)]
  if (length(va) > 100) va <- va[1:100]
  va <- scales::rescale(va, to = c(1, 10))
  vb <- lapply(names(va), function(x) names(igraph::neighbors(n, x)))
  vb <- unlist(unique(vb))
  va <- data.frame(V1 = names(va), value = va)
  va$x <- c(nrow(va):1)
  l <- list()
  for (i in vb) {
    z <- names(igraph::neighbors(n, i))
    z <- z[z %in% va$V1]
    if (length(z) < 1) next
    nnn <- va$value[va$V1 %in% z] / sum(va$value[va$V1 %in% z])
    l[[i]] <- data.frame(V1 = z, tot = nnn, V2 = i)
  }
  l <- dplyr::bind_rows(l)
  if (nrow(l) < 1) return(NULL)
  dplyr::distinct(merge(va, l, by = "V1"))
}
