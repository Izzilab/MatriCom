# On-screen plots and tables (from the original sub3(), sub4(), enplot(),
# infplot() and table observers); they receive the selection-filtered table.

blank_plot <- function() {
  err <- data.frame(x = 1, y = 1)
  pl <- ggplot(err, aes(x, y)) + geom_point()
  pl <- pl + theme_bw()
  pl <- pl + theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
                   axis.text.x = element_blank(), axis.text.y = element_blank(),
                   axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
                   panel.grid = element_blank(), plot.title = element_text(size = 20, face = "bold"))
  pl + labs(fill = '')
}

png_device <- function(filename, width, height, res) {
  if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png(filename, width = width, height = height, res = res)
  else grDevices::png(filename, width = width, height = height, res = res)
}

sub3_plot <- function(out) {
  d1 <- as.data.frame(table(out$`Type of interaction`))
  if (nrow(d1) < 1) return(ggplotly(ggplot()))
  d1$Freq <- round((d1$Freq/sum(d1$Freq))*100, 1)
  d1 <- merge(d1, data.frame(Var1 = c("cell-matrix", "extracellular")), by = "Var1", all.y = TRUE)
  d1$Freq[is.na(d1$Freq)] <- 0
  colnames(d1) <- c("type", "Percentage of total communications")
  d1$type <- ifelse(d1$type %in% "cell-matrix", "Cell-Matrisome", "Matrisome-Matrisome")
  d1$type <- ifelse(d1$type %in% "Cell-Matrisome", "Non.matrisome-Matrisome", d1$type)
  d1$type <- factor(d1$type, levels = c("Non.matrisome-Matrisome", "Matrisome-Matrisome"))
  pl <- ggplot(d1, aes(type, `Percentage of total communications`, fill = type)) +
    geom_bar(stat = "Identity") +
    scale_fill_manual(breaks = d1$type, values = c("#C5F385", "#83da07")) +
    ylim(c(0, 100))
  pl <- pl + theme_bw()
  pl <- pl + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
                   axis.ticks.x = element_blank(), panel.grid = element_blank(),
                   plot.title = element_text(size = 20, face = "bold"))
  pl <- pl + labs(fill = '') + ggtitle("Communication Pairs")
  ggplotly(pl, tooltip = c("Percentage of total communications"))
}

sub4_plot <- function(out, mlist) {
  out <- merge(out, mlist, by.x = "Gene1", by.y = "gene", all.x = T)
  out$family[is.na(out$family)] <- "Non.matrisome"
  out <- merge(out, mlist, by.x = "Gene2", by.y = "gene", all.x = T)
  out$family.y[is.na(out$family.y)] <- "Non.matrisome"
  df <- as.data.frame(table(out$family.x, out$family.y))
  df <- df[df$Freq > 0, ]
  if (nrow(df) < 1) return(ggplotly(ggplot()))
  lv <- c("ECM glycoproteins", "Collagens", "Proteoglycans", "ECM-affiliated proteins",
          "ECM regulators", "Secreted factors", "Non.matrisome")
  tb <- matrix(0, 7, 7, dimnames = list(lv, lv))
  for (i in 1:nrow(df)) {
    if (tb[rownames(tb) %in% df$Var2[i], colnames(tb) %in% df$Var1[i]] == 0) {
      tb[rownames(tb) %in% df$Var1[i], colnames(tb) %in% df$Var2[i]] <-
        sum(tb[rownames(tb) %in% df$Var1[i], colnames(tb) %in% df$Var2[i]], df$Freq[i])
    } else {
      tb[rownames(tb) %in% df$Var2[i], colnames(tb) %in% df$Var1[i]] <-
        sum(tb[rownames(tb) %in% df$Var2[i], colnames(tb) %in% df$Var1[i]], df$Freq[i])
    }
  }
  df <- reshape2::melt(tb)
  df <- df[df$value > 0, ]
  if (nrow(df) < 1) return(ggplotly(ggplot()))
  df <- data.frame(Var1 = as.character(df$Var1), Var2 = as.character(df$Var2), value = df$value)
  names(df)[3] <- "Freq"
  df$Freq <- as.numeric(df$Freq)
  df$Freq <- round((df$Freq/sum(df$Freq))*100, 3)
  names(df)[1:2] <- c("Gene1", "Gene2")
  df$Pair <- paste0(df$Gene1, "-", df$Gene2)
  names(df)[3] <- "Percentage of total communications"
  df <- distinct(merge(df, nd1, by.x = "Gene1", by.y = "V1", all.x = T))
  df$Gene1 <- factor(df$Gene1, levels = rev(lv))
  df$Gene2 <- factor(df$Gene2, levels = lv)
  pl <- ggplot(df, aes(Gene2, Gene1, label = Pair)) +
    geom_point(aes(size = `Percentage of total communications`, color = as.character(V2))) +
    scale_color_manual(breaks = df$V2, values = as.character(df$V2))
  pl <- pl + theme_bw()
  pl <- pl + theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
                   axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
                   panel.grid = element_blank(), plot.title = element_text(size = 20, face = "bold"))
  pl <- pl + labs(fill = '') + ggtitle("Matrisome Pairs") + theme(legend.position = "none")
  ggplotly(pl, tooltip = c("Pair", "Percentage of total communications"))
}

enrich_levels <- function(rr, s1) {
  s <- unique(rr$signature)
  s2 <- s[!(s %in% s1)]
  s2 <- s2[order(s2)]
  factor(rr$signature, levels = c(s1, s2))
}

enrich_plot <- function(rr) {
  rr$signature <- enrich_levels(rr, c("Core matrisome", "Matrisome-associated", "Collagens", "Proteoglycans",
                                      "ECM glycoproteins", "ECM regulators", "Secreted factors"))
  if (all(rr$p.value > 0.05)) return(ggplotly(ggplot()))
  rr <- rr[rr$p.value < 0.05, ]
  rr$p.value <- paste0("p value = ", round(rr$p.value, 4))
  colnames(rr) <- stringr::str_to_title(colnames(rr))
  gp <- ggplot(rr, aes(Signature, Populations, size = Overlap, color = Signature, text = P.value)) +
    geom_point() +
    scale_y_discrete(labels = abbreviate) +
    theme_bw() + xlab("") + ylab("") +
    theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(), plot.title = element_text(size = 20, face = "bold"),
          plot.title.position = "plot", legend.position = "none")
  gp <- gp + ggtitle("Matrisome-specific Signature Enrichment")
  ggplotly(gp, tooltip = c("Populations", "color", "Overlap", "P.value")) %>%
    layout(title = list(x = 0.01))
}

# table 3 / XLSX "enrichments" sheet; NULL when nothing is significant
enrich_table <- function(rr) {
  rr$signature <- enrich_levels(rr, c("Core matrisome", "Matrisome-associated", "Collagens", "Proteoglycans",
                                      "ECM glycoproteins", "ECM regulators", "Secreted factors"))
  if (all(rr$p.value > 0.05)) return(NULL)
  rr <- rr[rr$p.value < 0.05, ]
  rr$p.value <- paste0("p value = ", round(rr$p.value, 4))
  names(rr)[1:3] <- c("Populations", "Signature", "Overlap")
  rownames(rr) <- NULL
  rr
}

infl_plot <- function(va, mlist) {
  if (is.null(va)) return(ggplotly(ggplot()))
  names(va) <- c("Influencer", "size", "x", "Strength", "Influenced")
  va <- distinct(merge(va, mlist, by.x = "Influencer", by.y = "gene", all.x = T))
  va$category[is.na(va$category)] <- "Non.matrisome"
  va$color <- ifelse(va$category %in% "Core matrisome", "#002253",
                     ifelse(va$category %in% "Matrisome-associated", "#DB3E18", "grey80"))
  names(va)[6] <- "Influencer.Matrisome.Division"
  x <- unique(va$x)[order(unique(va$x))]
  names(x) <- c(1:length(x))
  df <- data.frame(x = x, new.x = names(x))
  va <- distinct(merge(va, df, by = "x"))
  va$new.x <- as.factor(va$new.x)
  gp <- ggplot(va, aes(new.x, Influenced, color = color, label = Influencer)) +
    geom_point(aes(size = Strength)) +
    scale_color_manual(breaks = as.character(va$color), values = as.character(va$color)) +
    theme_bw() + xlab("Influencers") + ylab("Influenced") +
    theme(axis.ticks = element_blank(), axis.text = element_blank(), panel.grid = element_blank(),
          plot.title = element_text(size = 20, face = "bold"), plot.title.position = "plot",
          legend.position = "none") +
    coord_flip()
  gp <- gp + ggtitle("Normalized Influence")
  ggplotly(gp, tooltip = c("Influencer", "Influenced", "Strength")) %>% layout(title = list(x = 0.01))
}

# main results table (on screen)
comm_table <- function(out, mlist) {
  out$cc <- NULL; out$g1g2 <- NULL; out$id <- NULL
  out <- distinct(merge(out, mlist, by.x = "Gene1", by.y = "gene", all.x = T))
  out <- distinct(merge(out, mlist, by.x = "Gene2", by.y = "gene", all.x = T))
  out[is.na(out)] <- "Non.matrisome"
  names(out)[11:14] <- c("Matrisome.Division.Gene1", "Matrisome.Category.Gene1",
                         "Matrisome.Division.Gene2", "Matrisome.Category.Gene2")
  out$perc.expr.Population1 <- round(out$perc.expr.Population1*100, 1)
  out$perc.expr.Population2 <- round(out$perc.expr.Population2*100, 1)
  out$mean.expr.Gene1 <- round(out$mean.expr.Gene1, 2)
  out$mean.expr.Gene2 <- round(out$mean.expr.Gene2, 2)
  out$Reliability.score <- as.numeric(out$Reliability.score)
  out <- out[order(-out$Reliability.score), ]
  out <- out[, c(2:6, 1, 7:ncol(out))]
  out <- out[, c(2, 1, 5, 6, 9, 10, 4, 8, 3, 7, 11:14)]
  out <- out[, c(1:10, 12, 14)]
  out$Matrisome.Category.Gene1 <- as.character(out$Matrisome.Category.Gene1)
  out$Matrisome.Category.Gene2 <- as.character(out$Matrisome.Category.Gene2)
  out$`Type of interaction` <- ifelse(out$`Type of interaction` == "cell-matrix", "Non.matrisome-Matrisome", "Matrisome-Matrisome")
  names(out)[3] <- "Type of communication"
  out
}

# "communication network" sheet of the XLSX export
comm_export <- function(out, mlist, CCgenes2) {
  out$cc <- NULL; out$g1g2 <- NULL
  out <- distinct(merge(out, mlist, by.x = "Gene1", by.y = "gene", all.x = T))
  out <- distinct(merge(out, mlist, by.x = "Gene2", by.y = "gene", all.x = T))
  out[is.na(out)] <- "Non.matrisome"
  names(out)[11:14] <- c("Matrisome.Division.Gene1", "Matrisome.Category.Gene1",
                         "Matrisome.Division.Gene2", "Matrisome.Category.Gene2")
  out$perc.expr.Population1 <- round(out$perc.expr.Population1*100, 1)
  out$perc.expr.Population2 <- round(out$perc.expr.Population2*100, 1)
  out$mean.expr.Gene1 <- round(out$mean.expr.Gene1, 2)
  out$mean.expr.Gene2 <- round(out$mean.expr.Gene2, 2)
  out$Reliability.score <- as.numeric(out$Reliability.score)
  out <- out[order(-out$Reliability.score), ]
  out <- out[, c(3, 2, 6, 1, 9, 10, 5, 8, 4, 7, 12, 11, 14, 13)]
  out[, 3] <- ifelse(out[, 3] %in% "cell-matrix", "Non.matrisome-Matrisome", "Matrisome-Matrisome")
  names(out)[3] <- "Type of communication"
  dd <- CCgenes2
  dd$label <- ifelse(dd$label == "extracellular", "extracellular (non-matrisome)", dd$label)
  out <- distinct(merge(out, dd, by.x = "Gene1", by.y = "gene"))
  out <- distinct(merge(out, dd, by.x = "Gene2", by.y = "gene"))
  names(out)[15:16] <- c("Compartment.Gene1", "Compartment.Gene2")
  out[, c(3, 2, 4, 1, 5, 6:12, 15, 13, 14, 16)]
}
