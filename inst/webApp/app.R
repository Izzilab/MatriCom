library(shiny)
library(Matrix)
library(DT)
library(dplyr)
library(data.table)
library(ggplot2)
library(shinyWidgets)
library(shinybusy)
library(shinyjs)
library(shinyhelper)
library(shinycssloaders)
library(shinyBS)
library(Seurat)
library(SeuratDisk)
library(SeuratObject)
library(scCustomize)
library(sp)
library(uwot)
library(shinydashboard)
library(scater)
library(googledrive)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(gridExtra)
library(writexl)
library(BiocManager)
library(packcircles)
library(ggiraph)
library(plotly)
library(viridis)
library(collapsibleTree)
library(igraph)
library(scales)
library(qs)
library(igraph)
library(stringr)
library(reshape2)
library(htmltools)
library(cicerone)

# Moved to R/, because stuff in `R` gets sourced automatically when in package/project
#source("./Data/get_expressed_genes_mod.R")

### GUIDED TOUR ###
guide <- Cicerone$
  new(allow_close = T)$
  step(
    "step1",
    "Upload your data",
    "To analyze your scRNA-seq data with MatriCom, use the “Browse” button to select and upload your RDS, QS, or H5AD file here."
  )$
  step(
    "step2",
    "Ensembl ID conversion",
    "If your dataset contains Ensembl IDs, you must activate this option to convert them to HGCN (human) or MGI (mouse) Gene Symbols."
  )$
  step(
    "step5",
    "Select population identifier column",
    "Choose the metadata column that contains your cell identity labels (e.g., cell type, cluster)."
  )$
  step(
    "step6",
    "Set minimum mean gene expression threshold",
    "Choose the minimum mean expression level a gene must have across the population to be considered in the analysis."
  )$
  step(
    "step7",
    "Set % positive population threshold",
    "Choose the minimum % of cells that must express a gene above the threshold you selected above for that gene to be considered in the analysis."
  )$
  step(
    "step10",
    "Run the analysis and download results",
    "When you are satisfied with your choices, run the analysis and download results!"
  )$
  step(
    "step9",
    "Post-run filters",
    "After analysis is complete, use these additional filters to refine the output you wish to display. Note that deselecting these filters may alter the biological relevance of your results. We strongly suggest keeping the default options active. See “?” buttons for more information."
  )
### ### ### ### ###

options(repos = BiocManager::repositories(),
        Seurat.object.assay.version = 'v3')

# This is processed by data-raw/data_prep.R from inst/extdata/ and saved in data/
# intlist <- readRDS("./Data/lr2.RDS")
# mlist <- readRDS("./Data/matrisome.list.RDS")[[1]]
# excl <- as.data.frame(fread("./Data/multimers.txt",sep="\t"))
# signs <- as.data.frame(fread("./Data/naba-1.csv"))
# CCgenes2 <- readRDS("./Data/CCgenes2.RDS")

nd1 <- data.frame(V1=c(
  "Drosophila matrisome",
  "Nematode-specific core matrisome",
  "Nematode-specific matrisome-associated",
  "Putative Matrisome",
  "Core matrisome",
  "Matrisome-associated",
  "Apical Matrix",
  "Cuticular Collagens",
  "Cuticlins",
  "ECM glycoproteins",
  "Collagens",
  "Proteoglycans",
  "ECM-affiliated proteins",
  "ECM regulators",
  "Secreted factors",
  "Non-matrisome"),
  V2=c("#B118DB",
       "#B118DB",
       "#741B47",
       "#B118DB",
       "#002253",
       "#DB3E18",
       "#B3C5FC",
       "#B3C5FC",
       "#BF9000",
       "#13349D",
       "#0584B7",
       "#59D8E6",
       "#F4651E",
       "#F9A287",
       "#FFE188",
       "#D9D9D9"
  )
)

matricom.slim <- function(obj,
                          group.column,
                          min.pct,
                          expr.filter,
                          lr,
                          mlist,
                          homo.policy=F,
                          conv=F
){

  #if(length(unique(obj@meta.data[,group.column]))==nrow(obj@meta.data)){stop()}

  lr <- lr
  mlist <- mlist
  # obj <- seurat.obj
  DefaultAssay(obj) <- "RNA"

  rownames(obj@assays$RNA@data) <- toupper(rownames(obj@assays$RNA@data)) #enforcing humanity
  rownames(obj@assays$RNA@counts) <- toupper(rownames(obj@assays$RNA@counts)) #enforcing humanity - just a backup

  if(length(colnames(obj@assays$RNA@data))<1){
    if(length(colnames(obj@assays$RNA@counts))<1){
      colnames(obj@assays$RNA@data) <- paste0("sample.",c(1:ncol(obj@assays$RNA@data)))
      colnames(obj@assays$RNA@counts) <- paste0("sample.",c(1:ncol(obj@assays$RNA@counts)))
      rownames(obj@meta.data) <- paste0("sample.",c(1:ncol(obj@assays$RNA@data)))
    }else{
      colnames(obj@assays$RNA@data) <- colnames(obj@assays$RNA@counts)
    }
  }

  if(length(rownames(obj@assays$RNA@data))<1){
    if(length(rownames(obj@assays$RNA@counts))<1){
      rownames(obj@assays$RNA@data) <- paste0("gene.",c(1:nrow(obj@assays$RNA@data)))
      rownames(obj@assays$RNA@counts) <- paste0("gene.",c(1:nrow(obj@assays$RNA@counts)))
    }else{
      rownames(obj@assays$RNA@data) <- rownames(obj@assays$RNA@counts)
    }
  }

  Idents(obj) <- as.factor(obj@meta.data[,group.column])
  #if(isTRUE(mode)){rel.mode <- "absolute"}else{rel.mode<-"relative"}

  if(isTRUE(conv)){
    gs <- ensembldb::select(EnsDb.Hsapiens.v86, keys=rownames(obj@assays$RNA@data), keytype = "GENEID", columns = c("SYMBOL","GENEID"), multivals="asNA")
    rn <- data.frame(GENEID=rownames(obj@assays$RNA@data),ord=c(1:length(rownames(obj@assays$RNA@data))))
    gs <- distinct(merge(rn,gs,by="GENEID",all.x=T))
    gs$SYMBOL[is.na(gs$SYMBOL)] <- paste0("not.annotated.gene_",sample(c(1:1000000),1))
    gs <- gs[order(gs$ord),]
    rownames(obj@assays$RNA@data) <- gs$SYMBOL

    gs <- ensembldb::select(EnsDb.Hsapiens.v86, keys=rownames(obj@assays$RNA@counts), keytype = "GENEID", columns = c("SYMBOL","GENEID"), multivals="asNA")
    rn <- data.frame(GENEID=rownames(obj@assays$RNA@counts),ord=c(1:length(rownames(obj@assays$RNA@counts))))
    gs <- distinct(merge(rn,gs,by="GENEID",all.x=T))
    gs$SYMBOL[is.na(gs$SYMBOL)] <- paste0("not.annotated.gene_",sample(c(1:1000000),1))
    gs <- gs[order(gs$ord),]
    rownames(obj@assays$RNA@counts) <- gs$SYMBOL
  }

  os <- obj@meta.data[,group.column][!is.na(obj@meta.data[,group.column])]

  obj <- subset(obj, features = unique(c(lr$V1,lr$V2)), idents=unique(os))
  obj@meta.data[,group.column] <- gsub("-","_",obj@meta.data[,group.column])
  obj@meta.data[,group.column] <- make.names(obj@meta.data[,group.column])
  Idents(obj) <- obj@meta.data[,group.column]

  obj <- subset(x = obj, idents = unique(na.omit(Idents(obj))))
  tab <- as.data.frame(table(obj@meta.data[,group.column]))
  #tab <- unique(tab[tab$Freq>min.pct,]$Var1)
  l <- nrow(tab)
  tab <- unique(tab[tab$Freq>1,]$Var1)

  if(length(tab)<l){
    obj <- subset(x = obj, idents = tab)
  }

  # remove unwanted assays
  # unw <- names(as.list(obj@assays))[names(as.list(obj@assays))!="RNA"]
  # if(length(unw)>=1){
  #   for(i in unw){
  #     obj[[i]] <- NULL
  #   }
  # }

  sender_celltypes <- unique(obj@meta.data[,group.column])
  list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes_mod, obj, (min.pct/100))

  #list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, obj, (1/100))

  names(list_expressed_genes_sender) <- sender_celltypes
  ligavals <- list()
  ligands <- list()
  for(i in names(list_expressed_genes_sender)){
    s <- obj@meta.data[,group.column]
    names(s) <- rownames(obj@meta.data)
    s <- names(s[s%in%i])
    z <- list_expressed_genes_sender[i][[1]]
    z <- z[z%in%unique(c(lr$V1,lr$V2))]
    k <- as.matrix(obj@assays$RNA@data[rownames(obj@assays$RNA@data)%in%z,
                                       colnames(obj@assays$RNA@data)%in%s])
    if(nrow(k)<2){next}else{
      if(ncol(k)<2){next}else{
        # k <- base::rowMeans(k)
        k <- rowMeans(replace(k, k == 0, NA), na.rm = TRUE)
        k <- k[k>=expr.filter]
        if(length(k)<1){next}else{
          ligavals[[i]] <- k
          ligands[[i]] <- names(k)
        }
      }}
  }

  ligands <- lapply(ligands,toupper)

  ints_omo <- list()
  for(i in names(ligands)){
    g <- ligands[[i]]
    omo <- lr[lr$V1%in%g & lr$V2%in%g,]
    if(nrow(omo)<1){next}else{
      omo$t1 <- ifelse(omo$V1%in%mlist$gene,"matrix","nonmatrix")
      omo$t2 <- ifelse(omo$V2%in%mlist$gene,"matrix","nonmatrix")
      omo$final <- paste0(omo$t1,"_",omo$t2)
      omo <- omo[omo$final!="nonmatrix_nonmatrix",]
      if(nrow(omo)<1){next}else{
        omo$t1 <- NULL
        omo$t2 <- NULL
        omo$final <- ifelse(omo$final=="matrix_matrix","extracellular","cell-matrix")
        omo$sender <- i
        omo$receiver <- i
        ints_omo[[i]] <- omo
      }
    }
  }
  ints_omo <- bind_rows(ints_omo)

  if(nrow(ints_omo)>0){
    ints_omo <- ints_omo[!duplicated(t(apply(ints_omo, 1, sort))),]
  }

  ints_etero <- list()
  for(i in names(ligands)){
    g <- ligands[[i]]
    lst <- list()
    for(w in names(ligands)[names(ligands)!=i])
      if(is.null(w)){next}else{
        g2 <- ligands[[w]]
        et_1 <- lr[lr$V1%in%g & lr$V2%in%g2,]
        et_2 <- lr[lr$V2%in%g & lr$V1%in%g2,]
        et <- bind_rows(et_1,et_2)
        if(nrow(et)<1){next}else{
          # et <- simplify(graph.data.frame(et))
          # et <- as.data.frame(as_edgelist(et))
          et$t1 <- ifelse(et$V1%in%mlist$gene,"matrix","nonmatrix")
          et$t2 <- ifelse(et$V2%in%mlist$gene,"matrix","nonmatrix")
          et$final <- paste0(et$t1,"_",et$t2)
          et <- et[et$final!="nonmatrix_nonmatrix",]
          if(nrow(et)<1){next}else{
            et$t1 <- NULL
            et$t2 <- NULL
            et$final <- ifelse(et$final=="matrix_matrix","extracellular","cell-matrix")
            et$sender <- i
            et$receiver <- w
            lst[[w]] <- et
          }
          ints_etero[[i]] <- bind_rows(lst)
        }
      }
  }

  ints_etero <- bind_rows(ints_etero)

  fin <- bind_rows(ints_omo,ints_etero)
  #fin[,3] <- NULL
  if(nrow(fin)<1){
    fin <- data.frame(t(rep(1,7)))
    # stop()
  }
  colnames(fin) <- c("Gene1","Gene2","source","relscore","Type of interaction","Population1","Population2")
  fin <- fin[,c(6,1,5,2,7,3,4)]


  if(isTruthy(homo.policy)){
    fin <- fin[!(fin$Population1==fin$Population2 & fin$Gene1==fin$Gene2), ]
  }else{
    fin <- fin
  }

  rownames(fin) <- NULL
  fin <- na.omit(fin)

  ### RELATIVE RELIABILITY SCORE (average % expressing cells per interaction pair * source reliability)

  # if(rel.mode=="relative"){

  if(nrow(fin)<1){
    fin <- data.frame(t(rep(1,7)))
    colnames(fin) <- c("Gene1","Gene2","source","relscore","Type of interaction","Population1","Population2")
    fin <- fin[,c(6,1,5,2,7,3,4)]
    fin$relscore <- 0
    names(fin)[7] <- "Reliability.score"
    return(fin)
  }else{
    # if(isTRUE(fin$Gene1=="1")){
    #   fin$source <- NULL
    #   fin$mean.expr.Gene1 <- 1
    #   fin$perc.expr.Population1 <- 1
    #   fin$mean.expr.Gene2 <- 1
    #   fin$perc.expr.Population2 <- 1
    #   fin <- fin[,c(1,2,7,8,3,4,9,10,5,6)]
    #   return(fin)
    # }else{
    prova <- Percent_Expressing(seurat_object = obj, features = unique(c(fin$Gene1,fin$Gene2)), group_by = group.column)
    prova <- prova/100
    # colnames(prova) <- gsub("\\."," ",colnames(prova))
    # fin$Population1 <- gsub("\\."," ",fin$Population1)
    # fin$Population2 <- gsub("\\."," ",fin$Population2)


    ndf <- apply(fin, 1, function(x){
      a <- prova[rownames(prova)%in%x[2],colnames(prova)%in%x[1]]
      b <- prova[rownames(prova)%in%x[4],colnames(prova)%in%x[5]]
      if(length(a)==0){a<-0}
      if(length(b)==0){b<-0}
      # v <- (mean(a,b) * as.numeric(x[7]))
      # return(v)
      # a <- round(a*100,1)
      # b <- round(b*100,1)
      x <- as.data.frame(t(x))
      x$perc.expr.Population1 <- a
      x$perc.expr.Population2 <- b
      x <- x[,c(1,2,8,3,4,9,5,6,7)]

      x <- x[x$perc.expr.Population1 >= (min.pct/100) & x$perc.expr.Population2 >= (min.pct/100), ]

      if(nrow(x)<1){
        x <- data.frame(t(rep("1",7)))
        colnames(x) <- c("Gene1","Gene2","source","relscore","Type of interaction","Population1","Population2")
        x[,4] <- "0"
        x <- x[,c(6,1,5,2,7,3,4)]
        #names(x)[7] <- "Reliability.score"
        x$perc.expr.Population1 <- 1
        x$perc.expr.Population2 <- 1
        x <- x[,c(1,2,8,3,4,9,5,6,7)]
        x$mean.expr.Gene1 <- 1
        x$mean.expr.Gene2 <- 1

        x <- x[,c(1,2,10,3,4,5,11,6,7:9)]

        return(x)

      }else{
        #names(x)[7] <- "Reliability.score"
        v1 <- ligavals[x$Population1][[1]]
        v1 <- v1[names(v1)%in%x$Gene1]
        v2 <- ligavals[x$Population2][[1]]
        v2 <- v2[names(v2)%in%x$Gene2]
        if(length(v1)==0){v1<-0}
        if(length(v2)==0){v2<-0}
        x$mean.expr.Gene1 <- v1
        x$mean.expr.Gene2 <- v2

        x <- x[,c(1,2,10,3,4,5,11,6,7:9)]

        # x <- x[x$mean.expr.Gene1 >= expr.filter & x$mean.expr.Gene2 >= expr.filter, ]

        return(x)
      }


    })
    # }

    ### ABSOLUTE RELIABILITY SCORE (average ranked expression per cell type per interaction pair * source reliability)

    # if(rel.mode=="absolute"){
    #   N <- nrow(obj@assays$RNA@data)
    #   M <- list()
    #   for(i in 1:ncol(obj@assays$RNA@data)){
    #     x <- obj@assays$RNA@data[,i]
    #     v <- x[order(-x)]
    #     nm <- names(v)
    #     v <- 1/c(1:N)
    #     df <- data.frame(gene=nm,value=v)
    #     M[[i]] <- df
    #   }
    #   M <- suppressWarnings(Reduce(function(x, y) merge(x, y, by="gene", all=TRUE), M))
    #   rownames(M) <- M$gene
    #   M$gene <- NULL
    #   names(M) <- colnames(obj@assays$RNA@data)
    #
    #   M <- as.data.frame(t(M))
    #   M$cell <- rownames(M)
    #   gdf <- data.frame(cell=rownames(obj@meta.data),type=obj@meta.data[,group.column])
    #   M <- distinct(merge(gdf,M,by="cell"))
    #   M$cell <- NULL
    #   M <- suppressWarnings(aggregate(M,by=list(M$type),mean))
    #   #M$Group.1 <- gsub("\\."," ",M$Group.1)
    #   M$type <- NULL
    #
    #   ndf <- apply(fin, 1, function(x){
    #     a <- M[M$Group.1%in%x[1],colnames(M)%in%x[2]]
    #     b <- M[M$Group.1%in%x[5],colnames(M)%in%x[4]]
    #     v <- (mean(a,b) * as.numeric(x[7]))
    #     return(v)
    #   })
    # }

    # ndf <- as.numeric(unlist(ndf))
    # ndf[is.na(ndf)] <- 0
    # fin$relscore <- ndf
    # fin$relscore <- round(fin$relscore,3)

    fin <- bind_rows(ndf)
    fin$relscore <- as.numeric(fin$relscore)
    fin <- fin[order(-fin$relscore),]
    names(fin)[11] <- "Reliability.score"

    fin$source <- NULL
    fin <- fin[fin$Reliability.score>0,]
    fin <- fin[fin$mean.expr.Gene1 >= expr.filter & fin$mean.expr.Gene2 >= expr.filter, ]
    fin <- fin[fin$perc.expr.Population1 >= (min.pct/100) & fin$perc.expr.Population2 >= (min.pct/100), ]
    # fin <- fin[fin$perc.expr.Population1 >= (input$cellprop/100) & fin$perc.expr.Population2 >= (input$cellprop/100), ]

    fin <- distinct(fin)

    return(fin)
    # }
  }


}

matrienrich <- function(data,fams,signs){
  # m <- data
  g <- fams[,c(2,1)]
  g2 <- fams[,c(3,1)]
  names(g) <- names(signs)
  names(g2) <- names(signs)
  l <- bind_rows(g,g2,signs)
  data$celltypes <- paste0(data$Population1,"-",data$Population2)

  enr <- list()
  for(i in unique(data$celltypes)){
    z <- data[data$celltypes==i,]
    group1 <- unique(c(z$Gene1,z$Gene2))
    ls <- list()
    for(w in unique(l$standard_name)){
      k <- unique(l[l$standard_name==w,2])
      o <- length(intersect(group1,k))
      p <- phyper(o-1, length(k), 4048-length(k), length(group1),lower.tail= FALSE)
      df <- data.frame(populations=i,signature=w,overlap=o,p.value=p)
      #if(df$p.value<0.05){df$overlap <- NA}
      ls[[w]] <- df
    }
    enr[[i]] <- bind_rows(ls)
  }
  enr <- bind_rows(enr)

  return(enr)
}

bbplot <- function(data){
  # m <- data

  data$comb <- paste0(data$Population1," - ",data$Population2)

  vertices <- as.data.frame(table(data$comb))
  names(vertices) <- c("name","size")
  #vertices$nnname <- c(1:nrow(vertices))
  vertices$shortName <- vertices$name

  packing <- circleProgressiveLayout(vertices$size, sizetype='area')
  v <- cbind(vertices, packing)
  v.gg <- circleLayoutVertices(packing, npoints=50)

  p <- ggplot() +
    geom_polygon_interactive(data = v.gg, aes(x, y, group = id, fill=id, tooltip = paste0(v$name[id]," : ",v$size[id]), data_id = v$name[id]), colour = "black", alpha = 0.6) +
    scale_fill_viridis() +
    geom_text(data = v, aes(x, y, label = size), size=2, color="black") +
    theme_void() +
    theme(legend.position="none",
          plot.title = element_text(size = 10, face = "bold"),
          #plot.margin=unit(c(0,0,0,0),"cm")
    ) +
    coord_equal() + ggtitle("Global Communication Cluster Map")

  css_default_hover <- girafe_css_bicolor(primary = "yellow", secondary = "red")

  return(suppressWarnings(girafe(ggobj = p,
                                 options = list(
                                   # opts_hover = opts_hover(css = css_default_hover),
                                   opts_hover_inv(css = "opacity:0.1;"),
                                   opts_hover(css = "stroke-width:2;",
                                              reactive = T),
                                   opts_zoom = opts_zoom(min = 1, max = 4),
                                   opts_tooltip = opts_tooltip(css = "padding:3px;background-color:#333333;color:white;"),
                                   opts_sizing = opts_sizing(rescale = TRUE),
                                   opts_toolbar = opts_toolbar(saveaspng = TRUE, position = "bottom", delay_mouseout = 1000),
                                   opts_selection(
                                     type = "single")
                                 ))))
} #for UI

savebbplot <- function(data){

  if(nrow(data)<1){
    err <- data.frame(x=1,y=1)
    pl <- ggplot(err,aes(x,y))
    pl <- pl + theme_bw()
    # pl <- pl + theme(legend.position = 'none')
    pl <- pl + theme(axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     axis.text.x = element_blank(),
                     axis.text.y = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.ticks.y = element_blank(),
                     panel.grid = element_blank(),
                     plot.title = element_text(size = 20, face = "bold"))
    # pl <- pl + scale_fill_viridis_d(option = "inferno")
    pl <- pl + labs(fill = '')
    return(pl)
  }else{

    data$comb <- paste0(data$Population1," - ",data$Population2)

    vertices <- as.data.frame(table(data$comb))
    names(vertices) <- c("name","size")
    vertices$nnname <- c(1:nrow(vertices))
    vertices$shortName <- paste0(vertices$name,"\n",vertices$size)

    packing <- circleProgressiveLayout(vertices$size, sizetype='area')
    v <- cbind(vertices, packing)
    v.gg <- circleLayoutVertices(packing, npoints=50)

    p <- ggplot() +
      geom_polygon_interactive(data = v.gg, aes(x, y, group = id, fill=id, tooltip = paste0(v$name[id]," : ",v$size[id]), data_id = v$name[id]), colour = "black", alpha = 0.6) +
      scale_fill_viridis() +
      geom_text(data = v, aes(x, y, label = shortName), size=2, color="black") +
      theme_void() +
      theme(legend.position="none",
            plot.title = element_text(size = 10, face = "bold"),
            #plot.margin=unit(c(0,0,0,0),"cm")
      ) +
      coord_equal() + ggtitle("Global Communication Cluster Map")

    return(p)
  }
}

savesub3 <- function(data,sel){
  out <- data
  out$cc <- paste0(out$Population1," - ",out$Population2)
  out <- out[out$cc %in% sel, ]
  if(nrow(out)<1){
    out <- data
    out$cc <- NULL
  }else{
    out$cc <- NULL
  }

  d1 <- as.data.frame(table(out$`Type of interaction`))
  if(nrow(d1)<1){
    err <- data.frame(x=1,y=1)
    pl <- ggplot(err,aes(x,y))
    pl <- pl + theme_bw()
    # pl <- pl + theme(legend.position = 'none')
    pl <- pl + theme(axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     axis.text.x = element_blank(),
                     axis.text.y = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.ticks.y = element_blank(),
                     panel.grid = element_blank(),
                     plot.title = element_text(size = 20, face = "bold"))
    # pl <- pl + scale_fill_viridis_d(option = "inferno")
    pl <- pl + labs(fill = '')
    return(pl)
  }else{
    d1$Freq <- round((d1$Freq/sum(d1$Freq))*100,1)
    nodf <- data.frame(Var1=c("cell-matrix","extracellular"))
    d1 <- merge(d1,nodf,all.y=T)
    d1$Freq[is.na(d1$Freq)] <- 0
    if(nrow(d1)<1){
      err <- data.frame(x=1,y=1)
      pl <- ggplot(err,aes(x,y)) + geom_point()
      pl <- pl + theme_bw()
      # pl <- pl + theme(legend.position = 'none')
      pl <- pl + theme(axis.title.x = element_blank(),
                       axis.title.y = element_blank(),
                       axis.text.x = element_blank(),
                       axis.text.y = element_blank(),
                       axis.ticks.x = element_blank(),
                       axis.ticks.y = element_blank(),
                       panel.grid = element_blank(),
                       plot.title = element_text(size = 20, face = "bold"))
      # pl <- pl + scale_fill_viridis_d(option = "inferno")
      pl <- pl + labs(fill = '')
      return(pl)
    }else{
      # d1$Freq.x[is.na(d1$Freq.x)] <- 0
      # d1$Freq.y <- NULL
      colnames(d1) <- c("type","Percentage of total communications")
      d1$type <- ifelse(d1$type%in%"cell-matrix","Cell-Matrisome","Matrisome-Matrisome")
      # if(nrow(d1)<2){
      #   if("cell-matrisome"%in%d1$type){
      #     d1 <- bind_rows(d1,
      #                     data.frame(type="matrisome-matrisome",`Percentage of total interactions`=0))
      #   }else{
      #     d1 <- bind_rows(d1,
      #                     data.frame(type="cell-matrisome",`Percentage of total interactions`=0))
      #   }
      # }
      d1$type <- factor(d1$type,levels=c("Cell-Matrisome","Matrisome-Matrisome"))

      pl <- ggplot(d1,aes(type,`Percentage of total communications`,fill=type)) +
        geom_bar(stat="Identity") +
        #scale_fill_nejm() +
        scale_fill_manual(breaks = d1$type, values = c("#C5F385","#83DA07"))+
        ylim(c(0,100)) #+
      #geom_vline(xintercept = 3, linetype=2, color="grey")


      pl <- pl + theme_bw()
      # pl <- pl + theme(legend.position = 'none')
      pl <- pl + theme(axis.title.x = element_blank(),
                       axis.text.x = element_blank(),
                       # axis.text.y = element_blank(),
                       axis.ticks.x = element_blank(),
                       panel.grid = element_blank(),
                       plot.title = element_text(size = 20, face = "bold"))
      # pl <- pl + scale_fill_viridis_d(option = "inferno")
      pl <- pl + labs(fill = '') + ggtitle("Communication Pairs")
      return(pl)
    }}
  # else{
  # return(ggplot())
  # }
}

savesub4 <- function(data,sel,mlist){
  out <- data
  out$cc <- paste0(out$Population1," - ",out$Population2)
  out <- out[out$cc %in% sel, ]
  if(nrow(out)<1){
    out <- data
    out$cc <- NULL
  }else{
    out$cc <- NULL
  }

  out <- merge(out,mlist,by.x="Gene1",by.y="gene",all.x=T)
  out$family[is.na(out$family)] <- "Non-matrisome"
  out <- merge(out,mlist,by.x="Gene2",by.y="gene",all.x=T)
  out$family.y[is.na(out$family.y)] <- "Non-matrisome"

  df <- as.data.frame(table(out$family.x,out$family.y))
  df <- df[df$Freq>0,]
  if(nrow(df)<1){
    err <- data.frame(x=1,y=1)
    pl <- ggplot(err,aes(x,y))
    pl <- pl + theme_bw()
    # pl <- pl + theme(legend.position = 'none')
    pl <- pl + theme(axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     axis.text.x = element_blank(),
                     axis.text.y = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.ticks.y = element_blank(),
                     panel.grid = element_blank(),
                     plot.title = element_text(size = 20, face = "bold"))
    # pl <- pl + scale_fill_viridis_d(option = "inferno")
    pl <- pl + labs(fill = '')
    return(pl)
  }else{
    tb <- matrix(0,7,7)
    rownames(tb) <- c("ECM glycoproteins",
                      "Collagens",
                      "Proteoglycans",
                      "ECM-affiliated proteins",
                      "ECM regulators",
                      "Secreted factors",
                      "Non-matrisome")

    colnames(tb) <- c("ECM glycoproteins",
                      "Collagens",
                      "Proteoglycans",
                      "ECM-affiliated proteins",
                      "ECM regulators",
                      "Secreted factors",
                      "Non-matrisome")

    for(i in 1:nrow(df)){
      if(tb[rownames(tb)%in%df$Var2[i],
            colnames(tb)%in%df$Var1[i]] == 0
      ){
        tb[rownames(tb)%in%df$Var1[i],
           colnames(tb)%in%df$Var2[i]] <- sum(
             tb[rownames(tb)%in%df$Var1[i],
                colnames(tb)%in%df$Var2[i]],
             df$Freq[i]
           )
      }else{
        tb[rownames(tb)%in%df$Var2[i],
           colnames(tb)%in%df$Var1[i]] <- sum(
             tb[rownames(tb)%in%df$Var2[i],
                colnames(tb)%in%df$Var1[i]],
             df$Freq[i]
           )
      }
    }

    df <- reshape2::melt(tb)
    df <- df[df$value>0,]

    if(nrow(df)<1){
      err <- data.frame(x=1,y=1)
      pl <- ggplot(err,aes(x,y))
      pl <- pl + theme_bw()
      # pl <- pl + theme(legend.position = 'none')
      pl <- pl + theme(axis.title.x = element_blank(),
                       axis.title.y = element_blank(),
                       axis.text.x = element_blank(),
                       axis.text.y = element_blank(),
                       axis.ticks.x = element_blank(),
                       axis.ticks.y = element_blank(),
                       panel.grid = element_blank(),
                       plot.title = element_text(size = 20, face = "bold"))
      # pl <- pl + scale_fill_viridis_d(option = "inferno")
      pl <- pl + labs(fill = '')
      return(pl)
    }else{
      df$variable <- NULL
      ddf <- apply(df,1,function(x){
        if(x[1]%in%"Not matrisome"){
          data.frame(Var1=x[2],Var2=x[1],value=x[3])
        }else{
          data.frame(Var1=x[1],Var2=x[2],value=x[3])
        }
      })
      df <- bind_rows(ddf)
      names(df)[3] <- "Freq"
      df$Freq <- as.numeric(df$Freq)
      df$Freq <- round((df$Freq/sum(df$Freq))*100,3)

      names(df)[1:2] <- c("Gene1","Gene2")
      df$Pair <- paste0(df$Gene1,"-",df$Gene2)
      names(df)[3] <- "Percentage of total communications"
      df <- distinct(merge(df,nd1,by.x="Gene1",by.y="V1",all.x=T))
      df$Gene1 <- factor(df$Gene1,
                         levels = rev(c("ECM glycoproteins",
                                        "Collagens",
                                        "Proteoglycans",
                                        "ECM-affiliated proteins",
                                        "ECM regulators",
                                        "Secreted factors",
                                        "Non-matrisome")))
      df$Gene2 <- factor(df$Gene2,
                         levels = c("ECM glycoproteins",
                                    "Collagens",
                                    "Proteoglycans",
                                    "ECM-affiliated proteins",
                                    "ECM regulators",
                                    "Secreted factors",
                                    "Non-matrisome"))
      pl <- ggplot(df,aes(Gene1,Gene2,label=Pair)) +
        geom_point(aes(size=`Percentage of total communications`,color=as.character(V2))) +
        scale_color_manual(breaks=df$V2,values=as.character(df$V2))
      #+
      # ylim(c(0,100)) +
      # scale_fill_manual(breaks = r$`Matrisome component`, values = as.character(r$V2))+
      # geom_vline(xintercept = 3, linetype=2, color="grey")

      pl <- pl + theme_bw() +xlab("Matrisome category of Gene1") + ylab("Matrisome category of Gene2")
      # pl <- pl + theme(legend.position = 'none')
      pl <- pl + theme(# axis.title.x = element_blank(),
        # axis.title.y = element_blank(),
        # axis.text.x = element_blank(),
        # axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 20, face = "bold"))
      # pl <- pl + scale_fill_viridis_d(option = "inferno")
      pl <- pl + labs(fill = '') + ggtitle("Matrisome Pairs") + NoLegend()
      return(pl)

    }
  }
}

saveploten <- function(data,sel,mlist,signs){
  out <- data
  out$cc <- paste0(out$Population1," - ",out$Population2)
  out <- out[out$cc %in% sel, ]
  if(nrow(out)<1){
    out <- data
    out$cc <- NULL
  }else{
    out$cc <- NULL
  }

  m <- out
  rr <- matrienrich(m,mlist,signs)
  #rr$alp <- ifelse(rr$p.value<0.05,0.8,0.1)
  s <- unique(rr$signature)
  s1 <- c("Core matrisome","Matrisome-associated",
          "Collagens","Proteoglycans","ECM glycoproteins",
          "ECM regulators","Secreted factors")
  s2 <- s[!(s %in% s1)]
  s2 <- s2[order(s2)]
  s <- c(s1,s2)
  rr$signature <- factor(rr$signature,levels = s)
  #rr$eval <- ifelse(rr$p.value<0.05,NA,rr$p.value)
  if(all(rr$p.value>0.05)){
    err <- data.frame(x=1,y=1)
    pl <- ggplot(err,aes(x,y))
    pl <- pl + theme_bw()
    # pl <- pl + theme(legend.position = 'none')
    pl <- pl + theme(axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     axis.text.x = element_blank(),
                     axis.text.y = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.ticks.y = element_blank(),
                     panel.grid = element_blank(),
                     plot.title = element_text(size = 20, face = "bold"))
    # pl <- pl + scale_fill_viridis_d(option = "inferno")
    pl <- pl + labs(fill = '')
    return(pl)
  }else{
    rr <- rr[rr$p.value<0.05,]
    rr$p.value <- paste0("p value = ",round(rr$p.value,4))

    gp <- ggplot(rr,aes(signature,populations,size=overlap,color=signature,text=p.value)) +
      geom_point() +
      scale_x_discrete(guide = guide_axis(angle = 90)) +
      #geom_point(aes(alpha=alp),show.legend = FALSE) +
      theme_bw() + xlab("Signatures") + ylab("Populations") +
      #guides(x =  guide_axis(angle = 90)) +
      theme(axis.title = element_blank(),
            # axis.text.x = element_blank(),
            #axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank(),
            plot.title = element_text(size = 20, face = "bold"),
            plot.title.position = "plot",
            legend.position = "none")
    gp <- gp + #labs(size = 'Overlap (all p<0.05)') +
      ggtitle("Matrisome-specific Signature Enrichment")

    return(gp)
  }
}

saveplotin <- function(data,sel,mlist){
  out <- data
  out$cc <- paste0(out$Population1," - ",out$Population2)
  out <- out[out$cc %in% sel, ]

  if(nrow(out)<1){
    out <- data
    out$cc <- NULL
  }else{
    out$cc <- NULL
  }

  # n <- distinct(data.frame(V1=out$Gene1,V2=out$Gene2))
  n <- distinct(out[,c(2,6)])
  if(nrow(n)<1){
    err <- data.frame(x=1,y=1)
    pl <- ggplot(err,aes(x,y))
    pl <- pl + theme_bw()
    # pl <- pl + theme(legend.position = 'none')
    pl <- pl + theme(axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     axis.text.x = element_blank(),
                     axis.text.y = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.ticks.y = element_blank(),
                     panel.grid = element_blank(),
                     plot.title = element_text(size = 20, face = "bold"))
    # pl <- pl + scale_fill_viridis_d(option = "inferno")
    pl <- pl + labs(fill = '')
    return(pl)
  }else{

    n <- simplify(graph.data.frame(n))

    #va <- page.rank(n)$vector
    va <- degree(n)
    va <- va[order(-va)]
    if(length(va)<=100){va<-va}else{va<-va[1:100]}
    va <- rescale(va,to=c(1,10))
    # vb <- names(neighbors(n,names(va)))
    vb <- lapply(names(va), function(x){
      names(neighbors(n,x))
    })
    vb <- unlist(unique(vb))

    va <- data.frame(V1=names(va),value=va)
    va$x <- c(nrow(va):1)
    l <- list()
    for(i in vb){
      z <- names(neighbors(n,i))
      z <- z[z%in%va$V1]
      if(length(z)<1){next}else{
        nnn <- va$value[va$V1%in%z]/sum(va$value[va$V1%in%z])
        l[[i]] <- data.frame(V1=z,tot=nnn,V2=i)
      }
    }
    l <- bind_rows(l)

    if(nrow(l)<1){
      err <- data.frame(x=1,y=1)
      pl <- ggplot(err,aes(x,y))
      pl <- pl + theme_bw()
      # pl <- pl + theme(legend.position = 'none')
      pl <- pl + theme(axis.title.x = element_blank(),
                       axis.title.y = element_blank(),
                       axis.text.x = element_blank(),
                       axis.text.y = element_blank(),
                       axis.ticks.x = element_blank(),
                       axis.ticks.y = element_blank(),
                       panel.grid = element_blank(),
                       plot.title = element_text(size = 20, face = "bold"))
      # pl <- pl + scale_fill_viridis_d(option = "inferno")
      pl <- pl + labs(fill = '')
      return(pl)
    }else{

      va <- distinct(merge(va,l,by="V1"))
      names(va) <- c("Influencer","size","x","Strength","Influenced")
      va <- distinct(merge(va,mlist,by.x="Influencer",by.y="gene",all.x=T))
      va$category[is.na(va$category)] <- "Non-matrisome"
      va$color <- ifelse(va$category%in%"Core matrisome","#002253",
                         ifelse(va$category%in%"Matrisome-associated","#DB3E18","grey80"))
      names(va)[6] <- "Influencer.Matrisome.Division"
      x <- unique(va$x)[order(unique(va$x))]
      names(x) <- c(1:length(x))
      df <- data.frame(x=x,new.x=names(x))
      va <- distinct(merge(va,df,by="x"))
      va$new.x <- as.factor(va$new.x)

      gp <- ggplot(va,aes(new.x,Influenced,color=color,
                          label=Influencer,
                          # label=Influencer.Matrisome.Division
      )) +
        geom_point(aes(size=Strength)) +
        scale_x_discrete(breaks = unique(va$new.x),labels = unique(va$Influencer)) +
        scale_y_discrete(guide = guide_axis(angle = 90)) +
        scale_color_manual(breaks = as.character(va$color), values = as.character(va$color)) +
        theme_bw() + xlab("Influencers") + ylab("Influenced") +
        theme(#axis.title = element_blank(),
          # axis.text.x = element_blank(),
          # axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(size = 20, face = "bold"),
          plot.title.position = "plot",
          legend.position = "none") +
        coord_flip()
      gp <- gp + #labs(size = 'Overlap (all p<0.05)') +
        ggtitle("Normalized Influence")
    }
  }
}

filtres <- function(data,
                    expr, #input$postsel
                    exc, #input$excl
                    type, #input$postsel2
                    hom, #input$nohomo
                    cc, #input$postsel3
                    dups #input$mmod
){

  m <- data


  if(isTruthy(all(m$Population1==1 & m$Population2==1 &nrow(m)==1))){
    return(m)
  }else{


    m <- m[m$Reliability.score%in%expr,]


    if(nrow(m)<1){
      df <- data.frame(error.message="No interactions found with these parameters. Please check your selection!")
      attributes(df)$outcome <- "failure"
      return(df)
    }

    if(isTruthy(hom)){
      m <- m[m$Gene1 != m$Gene2,]
      if(nrow(m)<1){
        df <- data.frame(error.message="No interactions found with these parameters. Please check your selection!")
        attributes(df)$outcome <- "failure"
        return(df)
      }else{
        m<-m
      }
    }

    if(isFALSE(exc)){
      m <- m
    }else{
      rskgn <- unique(unlist(strsplit(excl$multi,split=",")))
      rskgn <- gsub(" ","",rskgn)

      r1 <- m[m$Population1==m$Population2,] #safe

      r2 <- m[m$Population1!=m$Population2,] #potential risk
      r2b <- r2[!(r2$Gene1%in%rskgn) & !(r2$Gene2%in%rskgn),] #safe - to be added to r1 in the end
      r2a <- r2[r2$Gene1%in%rskgn & r2$Gene2%in%rskgn,] #very much risk - to be scored now!

      ex <- apply(excl,1,function(x){
        v <- unlist(strsplit(x,split=","))
        v <- gsub(" ","",v)
        return(v)
      })

      tst <- apply(r2a,1,function(x){
        a <- x[2]
        b <- x[4]
        l <- list()
        for(w in 1:length(ex)){
          z <- ex[[w]]
          l[[w]] <- ifelse(a%in%z & b%in%z,"remove","keep")
        }
        fn <- unlist(l)
        fn <- ifelse("remove"%in%fn,"remove","keep")
        return(fn)
      })

      r2a$final <- tst
      r2a <- r2a[r2a$final!="remove",]

      if(nrow(r2a)<1){
        clean <- bind_rows(r1,r2b)
        rownames(clean) <- NULL
        # return(clean)
      }else{
        clean <- bind_rows(r1,r2b,r2a)
        rownames(clean) <- NULL
        clean$final <- NULL
        # return(clean)
      }

      if(nrow(clean)<1){
        df <- data.frame(error.message="No interactions found with these parameters. Please check your selection!")
        attributes(df)$outcome <- "failure"
        return(df)
      }else{
        m <- clean
      }
    }

    if(length(type)>1){
      m <- m
    }else{
      if(type%in%"Homocellular"){
        m <- m[m$Population1==m$Population2,]
        if(nrow(m)<1){
          df <- data.frame(error.message="No interactions found with these parameters. Please check your selection!")
          attributes(df)$outcome <- "failure"
          return(df)
        }
      }else{
        m <- m[m$Population1!=m$Population2,]
        if(nrow(m)<1){
          df <- data.frame(error.message="No interactions found with these parameters. Please check your selection!")
          attributes(df)$outcome <- "failure"
          return(df)
        }
      }
    }

    if(isFALSE(dups)){
      m$g1g2 <- paste0(m$Gene1,"_",m$Gene2) #just to be sure it gets removed later on and doesn't mess up table structure
    }else{
      m$g1g2 <- paste0(m$Gene1,"_",m$Gene2)

      n <- as.data.frame(table(m$g1g2,m$Reliability.score))
      n <- n[n$Freq>0,]
      nn <- as.data.frame(table(n$Var1))
      ok <- nn$Var1[nn$Freq<2]

      m.1 <- m[m$g1g2%in%ok,]
      m.2 <- m[!(m$g1g2%in%ok),]
      l <- list()
      for(i in unique(m.2$g1g2)){
        z <- m.2[m.2$g1g2==i,]
        z <- z[order(-as.numeric(z$Reliability.score)),]
        l[[i]] <- z[1,]
      }
      m.2 <- bind_rows(l)
      m <- rbind(m.1,m.2)

      l2 <- list()
      for(i in unique(m$Reliability.score)){
        z <- m[m$Reliability.score==i,]
        z2 <- remove_reciprocalRows(z)
        if(nrow(z)<1){next}else{l2[[i]] <- z}
      }
      m <- bind_rows(l2)
    }

    CCgenes2$label <- str_to_sentence(CCgenes2$label)
    CCgenes2$label[CCgenes2$label=="Extracellular"] <- "Extracellular (Non-matrisome)"
    CCgenes2$label[CCgenes2$label=="Cell membrane"] <- "Intracellular"

    sg <- unique(CCgenes2[CCgenes2$label %in% cc,1])
    fin <- m[m$Gene1%in%sg & m$Gene2%in%sg,]
    if(nrow(fin)<1){
      df <- data.frame(error.message="No interactions found with these parameters. Please check your selection!")
      attributes(df)$outcome <- "failure"
      return(df)
    }else{
      attributes(fin)$outcome <- "success"
      return(fin)
    }
  }
}


# Define UI for application
ui <- fluidPage(

  use_cicerone(),

  tags$head(
    tags$style(HTML(
      "@keyframes glowing {
         0% { background-color: #fcfcfc; box-shadow: 0 0 5px #0795ab; }
         50% { background-color: #e8f0fc; box-shadow: 0 0 20px #43b0d1; }
         100% { background-color: #fcfcfc; box-shadow: 0 0 5px #0795ab; }
         }"
    ))),

    # Application title
    titlePanel("MatriCom"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
          useShinyjs(),
          useSweetAlert(),

          # Labs' ad!
          p("MatriCom developed by"),

          titlePanel(windowTitle = "MatriCom",

                     div(tags$a(img(
                       src = "izzilablogo.png",
                       width = "120px", height = "120px"),
                       # href="https://www.oulu.fi/en/research-groups/izzi-group",
                       # target="_blank"
                       ),tags$a(img(
                         src = "nabalablogo2.png",width = "240px", height = "120px"),
                         # href="https://sites.google.com/a/uic.edu/nabalab/",
                         # target="_blank"
                         )
                     )),
          p(em("If you use MatriCom in your publication, please cite our ",
               a("upcoming manuscript",
                 href = "https://pubmed.ncbi.nlm.nih.gov/",
                 target="_blank"))),

          p("Follow us at ",
            tags$img(
              src = "twitter.png",
              width = 20,
              height = 20
            ),
            a("the Matrisome Project",
              href = "https://twitter.com/Matrisome",
              target="_blank")
          ),

          a(actionButton(inputId = "email1", label = "Contact Us!",
                         icon = icon("envelope", lib = "font-awesome")),
            href="mailto:matrisomeproject@gmail.com"),

          actionButton("guide", "Take a guided tour of MatriCom",
                       icon = icon("circle-info", lib = "font-awesome",
                                   style = "animation: glowing 1300ms infinite;")),
          # Horizontal line
          tags$hr(),
          tags$h2(p("Data input")),
          tags$h3(p("Upload your own dataset")),
          div(
            id = "step1",
          # Input: Select a file
          fileInput("file1", "Select file",
                    multiple = F,
                    accept = c(".rds",
                               ".RDS",
                               ".H5AD",
                               ".h5ad",
                               ".qs",
                               ".QS"
                    )) %>% helper(type = "inline",
                                  title = "Accepted file types",
                                  content = c("MatriCom accepts scRNA-seq files of up to 1 GB from Seurat (RDS or QS format) and ScanPy/Loom (H5AD format).",
                                              " ",
                                              "IMPORTANT: MatriCom currently only accepts human and mouse datasets. If your file contains Ensembl gene IDs, you must convert them to HGCN (human) or MGI (mouse) Gene Symbols by activating the conversion button."),
                                  buttonLabel = "OK")),


          div(style = "margin-top: -17px"),

          div(
            id = "step2",
            awesomeCheckbox(
            inputId = "convbutton",
            label = "Convert my ENSEMBL IDs to gene symbols",
            value = FALSE,
            status = "danger"
          )),

          div(style = "margin-top: 30px"),

          tags$hr(),

          # Input: Select identities and other params
          tags$h2(p("Query Parameters")),
          div(
            id = "step5",
            uiOutput("id1")
            ),
          div(
            id = "step6",
            sliderInput("minexp", "Select the average gene expression value",
                      min = 1,
                      max = 10,
                      step = 1, value = 1) %>% helper(type = "inline",
                                                      title = "Mean gene expression threshold ",
                                                      content = c("This slider defines the minimum mean level above which each gene must be expressed, across the population, for that gene to be included in the analysis."),
                                                      buttonLabel = "OK")
             ),
          div(
            id = "step7",
          sliderInput("cellprop", "Select the minimum % positive population",
                        min = 1, max = 100, step = 1, value = 30) %>% helper(type = "inline",
                                                                             title = "Minimum % positive population",
                                                                             content = c("This slider defines the % of cells in a population that must be positive for a gene at the level chosen above for that gene to be included in the analysis."),
                                                                             buttonLabel = "OK")
          ),
          div(
            id = "step10",
            actionButton("runbutton", "run analysis", icon = icon("upload"),
                         style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),


            downloadButton("export",
                           label = "export all graphs (ZIP)",
                           class = "butt"),
            tags$head(tags$style(".butt{background:#cb9d06;} .butt{color: #f5f6f7;}")),

            downloadButton("export2",
                           label = "export tabular data (XLSX)",
                           class = "butt2"),
            tags$head(tags$style(".butt2{background:#74b72e;} .butt2{color: #f5f6f7;} .butt2:.noHover{
    pointer-events: none;}"))),
          tags$hr(),
          tags$h2(p("Filters")),
          div(
            id = "step9",


            materialSwitch(
              inputId = "mmod",
              label = "Maximize model",
              value = TRUE,
              status = "primary",
              right = TRUE) %>% helper(type = "inline",
                                       title = "Maximize model",
                                       content = c("Because the MatriCom omnibus ECM interaction database was curated from multiple sources, analysis may return duplicate entries of the same interaction with differing reliability scores. Selecting the ‘Maximize model’ option will only return interactions with the highest reliability score.",
                                                   " ",
                                                   "This filter also excludes \"reciprocal duplicates\", such that only one communication pair, defined here as the expression - by the same cell or two different cells - of two genes whose products can interact, is reported when the communicating populations and reliability scores are identical, i.e., if Population A expresses Gene A and Population B expresses Gene B, entries for Gene A-Gene B interactions will be treated the same as entries for Gene B-Gene A interactions."),
                                       buttonLabel = "OK"),

            materialSwitch(
              inputId = "excl",
              label = "Use exclusion list",
              value = TRUE,
              status = "primary",
              right = TRUE) %>% helper(type = "inline",
                                       title = "Exclude strictly homocellular interactions",
                                       content = c("Communications between some matrisome gene products are strictly homocellular, such that these multimeric proteins cannot be produced by the cooperation of multiple cell populations (e.g., collagen or laminin multimers must be assembled prior to secretion to the ECM).",
                                                   "When the ‘Use exclusion list’ option is selected, heterocellular multimers of this type are removed from the analysis results. We strongly suggest leaving this option active."),
                                       buttonLabel = "OK"),


            materialSwitch(
              inputId = "nohomo",
              label = "Remove homomeric interactions",
              value = TRUE,
              status = "primary",
              right = TRUE) %>% helper(type = "inline",
                                       title = "Exclude ambiguous homomeric interactions",
                                       content = c("Homocellular matrisome interactions involving multimers of the same protein (e.g., collagen I/collagen I in the same cell type) are difficult to score with scRNA-seq data.",
                                                   "When the ‘Remove homomeric interactions’ option is selected, all communication pairs of this type are removed from the analysis results. We strongly suggest leaving this option active."),
                                       buttonLabel = "OK"),

            prettyCheckboxGroup(
              inputId = "postsel",
              label = "Filter by reliability score",
              choices = c("3",
                          "2",
                          "1"),
              icon = icon("check"),
              status = "primary",
              inline = TRUE,
              selected = c("3"),
              animation = "jelly") %>% helper(type = "inline",
                                              title = "Reliability score",
                                              content = c("Users can filter results according to the reliability score of the returned communication pairs predicted to lead to protein interactions, with 1 being the lowest reliability and 3 being the highest.",
                                                          "Note that, by default, the initial output only displays interactions with a reliability score of 3."),
                                              buttonLabel = "OK"),

            prettyCheckboxGroup(
              inputId = "postsel2",
              label = "Filter by communication type",
              choices = c("Homocellular","Heterocellular"),
              icon = icon("check"),
              status = "primary",
              inline = TRUE,
              selected = c("Homocellular","Heterocellular"),
              animation = "jelly") %>% helper(type = "inline",
                                              title = "Communication type",
                                              content = c("Users can filter results according to the identity of cells expressing the communicating elements: homocellular (same population) or heterocellular (different populations)."),
                                              buttonLabel = "OK"),

            prettyCheckboxGroup(
              inputId = "postsel3",
              label = "Filter by cellular compartment",
              choices = c("Matrisome","Surfaceome","Extracellular (Non-matrisome)","Intracellular"),
              icon = icon("check"),
              status = "primary",
              inline = TRUE,
              selected = c("Matrisome","Surfaceome","Extracellular (Non-matrisome)"),
              animation = "jelly") %>% helper(type = "inline",
                                              title = "Cellular compartment",
                                              content = c("All communication pairs returned involve the product of at least one matrisome gene. Users may filter results according to the localization of the communicating protein: matrisome, cell surface (surfaceome), extracellular space (non-matrisome), or intracellular space.",
                                                          "In case a protein has multiple locations, we implemented the following hierarchy: Matrisome > Surfaceome > Extracellular (non-matrisome) > Intracellular. Genes encoding proteins that do not fall into one of the three highest-ranked compartments are marked as intracellular, which is deselected by default.",
                                                          " ",
                                                          "Sources for localization annotations are The Matrisome Project [DOI: 10.1074/mcp.M111.014647], the in silico human Surfaceome [DOI: 10.1073.pnas1808790115], and Gene Ontology [GO:0005576 (extracellular region)]."),
                                              buttonLabel = "OK")


          ),


        ), #END OF PANEL

        mainPanel(

          tabsetPanel(
            id = "tabset",
            type = "tabs",

            tabPanel("COMMUNICATION NETWORK",
                     div(style = "height: 20px"),
          fluidRow(
            column(width = 12, withSpinner(girafeOutput("plot_gir1",height = 1000
                                                        ),type=1)),

          ),

          div(style = "height: 50px"),
          fluidRow(
            column(width = 6, withSpinner(plotlyOutput("plot_sub3",height = 350),type=1)),
            column(width = 6, withSpinner(plotlyOutput("plot_sub4",height = 350),type=1)),
          ),
          div(style = "height: 50px"),
          DTOutput('tbl')
            ),

          tabPanel("NETWORK INFLUENCERS",
                   div(style = "height: 20px"),
          # div(style = "height: 50px"),
          fluidRow(
            column(width = 12,plotlyOutput("plot_in",height = 350),type=1),
            div(style = "height: 50px"),
            DTOutput('tbl2')
          ),


          ),

          tabPanel("ENRICHMENT ANALYSIS",
                   div(style = "height: 20px"),
          fluidRow(
            column(width = 12, plotlyOutput("plot_en",height = 1000),type=1),
            div(style = "height: 50px"),
            DTOutput('tbl3')
          ),
          ),
          
          tabPanel("ABOUT MATRICOM",
                   value = "help", 
                   # icon = tags$img(src='icons/about.svg',  
                   # height='48', width='48'
                   htmltools::includeMarkdown("www/about.md")
          ),
          # htmltools::includeMarkdown(paste0("www", "/","about.md"))
          
          
          )
        )
    )
)

# Define server logic ----

options(shiny.maxRequestSize=1000000*1024^2)

server <- function(input,output,session) {



  observe_helpers()
  guide$init()


  observeEvent(input$guide, {
    guide$start()
  }) #initiates the guide

  selected_state <- reactive({
    input$plot_gir1_selected
  }) #reactive to capture the selection from the bubble plot


  observeEvent(selected_state(),{
    if(!is.null(selected_state())){
      showNotification(
        paste0("You are now viewing results for ",selected_state(),".\nYou can remove the selection by clicking it again!"),
        duration = 4,
        closeButton = T,
        type = "message")
    }

  }) #side notes for bubble selection

  observeEvent(du(), {
    #req(du())
    nm <- sample(rownames(du()@assays$RNA@data),10)
      nm <- nm[!(grepl("orf",nm,ignore.case = T))]
      if(all(nm==toupper(nm))=="FALSE"){
        req(du())
        if(value2()==56){value2(1)}else{value2(56)}
        value3(0)
        showModal(modalDialog(
                title = "IMPORTANT!",
                "MatriCom has automatically detected that your input is likely not human. Please note that, in this case, MatriCom directly converts all genes to uppercase before proceeding with the analysis. This should generally work, but if you want to minimize the risk of loosing some genes, please convert your gene IDs before resubmitting your data!",
                easyClose = FALSE,
                footer = modalButton("close and continue")
              ))
      }else{
        req(!is.null(input$file1$datapath))
          show_modal_spinner(
            spin = "self-building-square",
            color = "grey",
            text = "Getting your file ready..."
          )
          req(du())
          if(value2()==56){value2(1)}else{value2(56)}
          value3(0)
          remove_modal_spinner()
      }
  }) #user-specific file upload spinner ### ALSO CONTAINS THE UPDATER TO VALUE2()!!! ###ALSO WARNS ABOUT MOUSE DATA!!!


  output$id1 <- renderUI({
    selectizeInput("ids", "Cell identity labels",
                   choices = "",
                   options = list(create = F)) %>% helper(type = "inline",
                                                          title = "select cell identities",
                                                          content = c("Datasets should contain at least one metadata column that holds cell identity labels (e.g., cell type, cluster, etc.) that will be used to define cell populations for analysis.",
                                                                      "Once a file is successfully uploaded, this menu will automatically populate with a list of the column headers."),
                                                          buttonLabel = "OK")
  }) #cell ID button


  #disable all buttons at start - they will pop up when a file is loaded
  shinyjs::disable("minexp")
  shinyjs::disable("cellprop")
  shinyjs::disable("scbutton")
  shinyjs::disable("mmod")
  shinyjs::disable("excl")
  shinyjs::disable("nohomo")
  shinyjs::disable("postsel")
  shinyjs::disable("postsel2")
  shinyjs::disable("postsel3")
  shinyjs::disable("runbutton")
  shinyjs::disable("export")
  shinyjs::disable("export2")

  value <- reactiveVal(0) #reactive value to store changes in data. Needed to switch samples

  value2 <- reactiveVal(56) #reactive value to keep track of user/OA switching and govern buttons accordingly

  value3 <- reactiveVal(0) #reactive value to control only the download buttons, especially if data are changed!

  du <- reactive({
    req(input$file1)

    dp <- input$file1$datapath

    if(isTruthy(length(dp[grepl(".RDS",dp,ignore.case = T)])>0)){
      df <- readRDS(input$file1$datapath)

      if(strsplit(as.character(df@version),split="\\.")[[1]][1] != 3){
        def <- DefaultAssay(df)
        mat <- df@assays[[def]]$counts
        rownames(mat) <- rownames(df@assays[[def]]$counts)
        colnames(mat) <- colnames(df@assays[[def]]$counts)
        mat <- CreateAssayObject(counts=mat)
        mat2 <- CreateSeuratObject(counts=mat,meta.data = df@meta.data)
        mat2 <- NormalizeData(mat2)
        #df[["RNA"]] <- as(object = df[["RNA"]], Class = "Assay")
      }

      value(1)
      return(mat2)
    }

    if(isTruthy(length(dp[grepl(".qs",dp,ignore.case = T)])>0)){
      df <- qread(input$file1$datapath)

      if(strsplit(as.character(df@version),split="\\.")[[1]][1] != 3){
        def <- DefaultAssay(df)
        mat <- df@assays[[def]]$counts
        rownames(mat) <- rownames(df@assays[[def]]$counts)
        colnames(mat) <- colnames(df@assays[[def]]$counts)
        mat <- CreateAssayObject(counts=mat)
        mat2 <- CreateSeuratObject(counts=mat,meta.data = df@meta.data)
        mat2 <- NormalizeData(mat2)
        #df[["RNA"]] <- as(object = df[["RNA"]], Class = "Assay")
      }

      value(1)
      return(mat2)
    }

    if(isTruthy(length(dp[grepl(".H5AD",dp,ignore.case = T)])>0)){
      #g <- gsub("(.*/\\s*(.*$))", "\\2", dp)
      g <- dp
      Convert(g, dest = "h5seurat", overwrite = TRUE, verbose = F)
      g <- gsub(".h5ad","",g)
      g <- gsub(".H5AD","",g)
      df <- LoadH5Seurat(paste0(g,".h5seurat"),verbose = F)

      if(strsplit(as.character(df@version),split="\\.")[[1]][1] != 3){
        def <- DefaultAssay(df)
        mat <- df@assays[[def]]$counts
        rownames(mat) <- rownames(df@assays[[def]]$counts)
        colnames(mat) <- colnames(df@assays[[def]]$counts)
        mat <- CreateAssayObject(counts=mat)
        mat2 <- CreateSeuratObject(counts=mat,meta.data = df@meta.data)
        mat2 <- NormalizeData(mat2)
        #df[["RNA"]] <- as(object = df[["RNA"]], Class = "Assay")
      }

      value(1)
      return(mat2)
    }


  }) #user-specific FILE READ-IN

  observeEvent(input$file1, {
    shinyjs::enable("runbutton")
  }) #activates the run button at file upload

  observeEvent(input$imprt, {
    shinyjs::enable("runbutton")
  }) #activates the run button at OA file import

  d0 <- reactive({
    if(value()==1){
      #req(du())
      d0 <- du()
    }else{
      if(value()==2){
        d0 <- dats()
      }else{
        # d0 <- NULL
        value()==0
      }
    }
    return(d0)
    }) #SWITCH BETWEEN USER UPLOAD AND OPEN-ACCESS

  observeEvent(du(), {
    nm1 <- names(du()@assays)
    if(length(intersect(nm1,"RNA"))<1){
      value(0)
      # value3(0)
      showNotification(
        "File does not have the expected \"RNA\" slot! You cannot proceed with analysis!",
        type = "error")
    }else{
      if(nrow(du()@assays$RNA@data)<1){
        value(0)
        # value3(0)
        showNotification(
          "File does not have the expected normalized counts! You cannot proceed with analysis!",
          type = "error")
      }
    }
  }) #side notes for file upload failure/wrong format

  observeEvent(du(),{
    if(value()==1){
         showNotification(
          "File correctly read! You can now proceed with analysis!",
          type = "message")
    }
  }) #side notes for file upload success/right format

  observeEvent(input$imprt, {
    req(dats())
    value(2)
    if(value2()==56){value2(1)}else{value2(56)}
    # value3(0)
    showNotification(
      "File correctly read! You can now proceed with analysis!",
      type = "message")
  }) #side notes for open access file upload success

  observeEvent(input$runbutton,{
    if(value2()==56){value2(1)}else{value2(56)}
    # value3(1)
  }) #value2 controller

  observe({
    req(!is.null(d0()))
    req(value()!=0 | value2()!=56)
    # if(value()!=0 & value2()==56){
      # shinyjs::enable("mmod")
      shinyjs::enable("minexp")
      shinyjs::enable("cellprop")
      shinyjs::enable("scbutton")
      # shinyjs::enable("runbutton")
      # shinyjs::disable("export")
      # shinyjs::disable("export2")
    # }
  }) #show buttons when the file is correct or changed

  observeEvent(input$file1, {
    shinyjs::disable("export")
    shinyjs::disable("export2")
  }) #deactivates the download buttons at file upload

  observeEvent(input$imprt, {
    shinyjs::disable("export")
    shinyjs::disable("export2")
  }) #deactivates the download buttons at file import

  observeEvent(input$runbutton, {
    req(d4())
    shinyjs::enable("export")
    shinyjs::enable("export2")
  }) #activates the download buttons at run

  # observe({
  # req(value3()!=0)
  #   shinyjs::enable("export")
  #   shinyjs::enable("export2")
  # }) #show download buttons only when the input is correct and a run has been made

  observe({
      req(value()!=0)
      req(d0())
      if(value()==1){
        nmv <- colnames(d0()@meta.data)
      }
      if(value()==2){
        #"Tabula Sapiens","The Human Protein Atlas","HuBMAP (Azimuth demo)"
        # if(input$crt=="Tabula Sapiens"){
        #   nmv <- "cell_type"
        #   #,"cell_type_ontology_term_id","disease","development_stage","sex","ethnicity")
        # }
        # if(input$crt=="The Human Protein Atlas"){
        #   nmv <- "cell_type"
        # }
        # if(input$crt=="HuBMAP"){
        #   nmv <- "cell_type"
        #   #nmv <- colnames(d0()@meta.data)
        # }
        nmv <- "cell_type"
       }

      updateSelectizeInput(session, "ids", choices = nmv, server = TRUE,
                           options = list(maxOptions=2000))
  }) #group column name button updater

  observeEvent(input$runbutton, {
    show_modal_spinner(
      spin = "trinity-rings",
      color = "blue",
      text = "Analysis has started..."
    )
    req(d4())
    value2(1)
    remove_modal_spinner()
  })#spinner for the initial wait... could be useful if the server is under load

  observe({
    req(d4())
    req(attributes(d2())$outcome != "failure")
    if(value()!=0 & value2()!=56){
      shinyjs::enable("excl")
      shinyjs::enable("nohomo")
      shinyjs::enable("postsel")
      shinyjs::enable("postsel2")
      shinyjs::enable("postsel3")
      shinyjs::enable("mmod")
    }else{
      shinyjs::disable("mmod")
      shinyjs::disable("excl")
      shinyjs::disable("nohomo")
      shinyjs::disable("postsel")
      shinyjs::disable("postsel2")
      shinyjs::disable("postsel3")
    }
  }) #show the filter buttons when interactions are found

  # observe({
  #   req(d4())
  #   shinyjs::enable("export")
  #   shinyjs::enable("export2")
  # }) #enable the download buttons when the final table is ready

  observe({
    req(d0())
    if(value2()==56){
      shinyjs::hide("plot_gir1")
      shinyjs::hide("plot_in")
      shinyjs::hide("plot_in2")
      shinyjs::hide("tree")
      shinyjs::hide("plot_sub3")
      shinyjs::hide("plot_sub4")
      shinyjs::hide("plot_en")
      shinyjs::hide("tbl")
      shinyjs::hide("tbl2")
      shinyjs::hide("tbl3")
      }else{
      shinyjs::show("plot_gir1")
      shinyjs::show("plot_in")
      shinyjs::show("plot_in2")
      shinyjs::show("tree")
      shinyjs::show("plot_sub3")
      shinyjs::show("plot_sub4")
      shinyjs::show("plot_en")
      shinyjs::show("tbl")
      shinyjs::show("tbl2")
      shinyjs::show("tbl3")
    }
  }) #hide graphs if value2 is zeroed

  observeEvent(input$postsel2,{
    req(d4())
    if(isTruthy(length(input$postsel2)<1)){
      showNotification(
        "WARNING: Deselecting both the communication types removes all interactions!",
        type = "warning")
    }
  }) #side notes for warning on filter out all interaction types

  observeEvent(input$postsel3,{
    req(d4())
    if(length(intersect("Matrisome",input$postsel3))<1){
      showNotification(
        "WARNING: Deselecting \"Matrisome\" removes all interactions!",
        type = "error")
    }
  }) #side notes for warning on filter out matrisome


  d2 <- eventReactive(input$runbutton,{
  req(d0())

    if(value()==1){

    res <- tryCatch({
              df <- matricom.slim(d0(),
                            input$ids,
                            input$cellprop,
                            input$minexp,
                            intlist,
                            mlist,
                            input$nohomo,
                            conv=F)

        # if(df$Reliability.score==0){
        if(nrow(df)<1){
          df <- data.frame(error.message="No interactions found with these parameters. Please check your selection!")
          attributes(df)$outcome <- "failure"
          return(df)
        }else{
          attributes(df)$outcome <- "success"
          return(df)
        }

    }, error = function(e) {
      #debug_msg(e$message)
      df <- data.frame(error.message="No interactions found with these parameters. Please check your selection!")
      attributes(df)$outcome <- "failure"
      return(df)
    })
    }

    if(value()==1 & isTruthy(input$convbutton)){

      res <- tryCatch({
        df <- matricom.slim(d0(),
                            input$ids,
                            input$cellprop,
                            input$minexp,
                            intlist,
                            mlist,
                            input$nohomo,
                            conv=TRUE)
        if(nrow(df)<1){
          df <- data.frame(error.message="No interactions found with these parameters. Please check your selection!")
          attributes(df)$outcome <- "failure"
          return(df)
        }else{
          attributes(df)$outcome <- "success"
          return(df)
        }
      }, error = function(e) {
        #debug_msg(e$message)
        df <- data.frame(error.message="No interactions found with these parameters. Please check your selection!")
        attributes(df)$outcome <- "failure"
        return(df)
      })
    }


    if(value()==2){
      res <- tryCatch({
        df <- d0()
        df <- df[df$perc.expr.Population1 >= (input$cellprop/100) & df$perc.expr.Population2 >= (input$cellprop/100) & df$mean.expr.Gene1>=input$minexp & df$mean.expr.Gene2>=input$minexp, ]


        # if(df$Reliability.score==0){
        if(nrow(df)<1){
          df <- data.frame(error.message="No interactions found with these parameters. Please check your selection!")
          attributes(df)$outcome <- "failure"
          return(df)
        }else{
          attributes(df)$outcome <- "success"
          return(df)
        }

      }, error = function(e) {
        #debug_msg(e$message)
        df <- data.frame(error.message="No interactions found with these parameters. Please check your selection!")
        attributes(df)$outcome <- "failure"
        return(df)
      })
    }

  }) #actual trycatch is here

  d4 <- reactive({
    req(d2())
    if(attributes(d2())$outcome != "failure"){
      m <- filtres(d2(),input$postsel,input$excl,input$postsel2,input$nohomo,input$postsel3,dups = input$mmod)
      return(m)
      #m <- m[m$perc.expr.Population1 >= (input$cellprop/100) & m$perc.expr.Population2 >= (input$cellprop/100) & m$mean.expr.Gene1>=input$minexp & m$mean.expr.Gene2>=input$minexp, ]
    }else{
      return(d2())
    }
  }) #filtering results according to postrun


  observe(
    if(attributes(d2())$outcome != "failure"){

      out <- d4()

      if(nrow(out)<1 | attributes(out)$outcome=="failure"){
        output$tbl <- renderDT(data.frame(error.message="No interactions found with these parameters. Please check your selection!"))
      }else{

      out$cc <- paste0(out$Population1," - ",out$Population2)
      out <- out[out$cc %in% selected_state(), ]
      if(nrow(out)<1){
        out <- d4()
        out$cc <- NULL
        out$g1g2 <- NULL
        out$id <- NULL
        # out <- out
        out <- distinct(merge(out,mlist,by.x="Gene1",by.y="gene",all.x=T))
        out <- distinct(merge(out,mlist,by.x="Gene2",by.y="gene",all.x=T))
        out[is.na(out)] <- "Non.matrisome"
        names(out)[11:14] <- c("Matrisome.Division.Gene1","Matrisome.Category.Gene1","Matrisome.Division.Gene2","Matrisome.Category.Gene2")
        # out$Reliability.score <- as.numeric(out$Reliability.score)
        # out <- out[order(-out$Reliability.score),]

      }else{
        out$cc <- NULL
        out$g1g2 <- NULL
        out$id <- NULL
        out <- distinct(merge(out,mlist,by.x="Gene1",by.y="gene",all.x=T))
        out <- distinct(merge(out,mlist,by.x="Gene2",by.y="gene",all.x=T))
        out[is.na(out)] <- "Non.matrisome"
        names(out)[11:14] <- c("Matrisome.Division.Gene1","Matrisome.Category.Gene1","Matrisome.Division.Gene2","Matrisome.Category.Gene2")

      }

      out$perc.expr.Population1 <- round(out$perc.expr.Population1*100,1)
      out$perc.expr.Population2 <- round(out$perc.expr.Population2*100,1)
      out$mean.expr.Gene1 <- round(out$mean.expr.Gene1,2)
      out$mean.expr.Gene2 <- round(out$mean.expr.Gene2,2)
      out$Reliability.score <- as.numeric(out$Reliability.score)
      out <- out[order(-out$Reliability.score),]
      out <- out[,c(2:6,1,7:ncol(out))]

      #change to table struct
      out <- out[,c(2,1,5,6,9,10,4,8,3,7,11:14)]
      # out[,12] <- out[,11]
      out <- out[,c(1:10,12,14)]
      out$Matrisome.Category.Gene1 <- as.character(out$Matrisome.Category.Gene1)
      out$Matrisome.Category.Gene2 <- as.character(out$Matrisome.Category.Gene2)
      out$`Type of interaction` <- ifelse(out$`Type of interaction`=="cell-matrix","Non.matrisome-Matrisome","Matrisome-Matrisome")
      names(out)[3] <- "Type of communication"

      output$tbl <- renderDT(DT::datatable(out, extensions = 'Buttons',
                             options = list(#scrollX=TRUE, #lengthMenu = c(5,10,15),
                                            #paging = TRUE, searching = TRUE,
                                            #fixedColumns = TRUE, autoWidth = TRUE,
                                            #ordering = TRUE, dom = 'tB',
                                            buttons = c('csv', 'excel','pdf')),
                             rownames = F)  %>%
                               formatStyle('Reliability.score',
                                           background = styleColorBar(c(0,3), 'lightblue'),
                                           backgroundSize = '95% 80%',
                                           backgroundRepeat = 'no-repeat',
                                           backgroundPosition = 'left')
                             %>%
                               formatStyle('perc.expr.Population1',
                                           background = styleColorBar(c(0,100), 'lightcoral'),
                                           backgroundSize = '95% 80%',
                                           backgroundRepeat = 'no-repeat',
                                           backgroundPosition = 'left')
                             %>%
                               formatStyle('perc.expr.Population2',
                                           background = styleColorBar(c(0,100), 'lightcoral'),
                                           backgroundSize = '95% 80%',
                                           backgroundRepeat = 'no-repeat',
                                           backgroundPosition = 'left')
                             %>%
                               formatStyle('mean.expr.Gene1',
                                           background = styleColorBar(c(0,max(out$mean.expr.Gene1)), 'bisque'),
                                           backgroundSize = '95% 80%',
                                           backgroundRepeat = 'no-repeat',
                                           backgroundPosition = 'left')
                             %>%
                               formatStyle('mean.expr.Gene2',
                                           background = styleColorBar(c(0,max(out$mean.expr.Gene2)), 'bisque'),
                                           backgroundSize = '95% 80%',
                                           backgroundRepeat = 'no-repeat',
                                           backgroundPosition = 'left'))


      }
    }else{
      output$tbl <- renderDT(data.frame(error.message="No interactions found with these parameters. Please check your selection!")
      )
    }
  ) #table

  # output$plot_gir1 <- renderEcharts4r({
  output$plot_gir1 <- renderGirafe({
    req(d4())
    req(value2()!=0)
    # withProgress(message = 'generating the interaction map', value = 0, {
    if(attributes(d4())$outcome != "failure"){
      out <- d4()
      if(nrow(out)<1){
        return(ggplot())
      }else{
        # bbplot2(d4())
        bbplot(d4())
      }
    }else{
      return(ggplot())
    }
    # incProgress(1)})
  }) #render the bubbleplot 2

  sub3 <- function(){
    req(d4())
    req(value2()!=0)
    if(attributes(d4())$outcome != "failure"){

      out <- d4()
      out$cc <- paste0(out$Population1," - ",out$Population2)
      out <- out[out$cc %in% selected_state(), ]
      if(nrow(out)<1){
        out <- d4()
        out$cc <- NULL
      }else{
        out$cc <- NULL
      }

      d1 <- as.data.frame(table(out$`Type of interaction`))
      if(nrow(d1)<1){
        return(ggplot())
      }else{
        d1$Freq <- round((d1$Freq/sum(d1$Freq))*100,1)
        nodf <- data.frame(Var1=c("cell-matrix","extracellular"),Freq=c(0,0))
        d1 <- merge(d1,nodf,by="Var1")
        d1$Freq.x[is.na(d1$Freq.x)] <- 0
        d1$Freq.y <- NULL
        colnames(d1) <- c("type","Percentage of total communications")
        d1$type <- ifelse(d1$type%in%"cell-matrix","Cell-Matrisome","Matrisome-Matrisome")
        if(nrow(d1)<2){
          if("cell-matrisome"%in%d1$type){
            d1 <- bind_rows(d1,
                            data.frame(type="Matrisome-Matrisome",`Percentage of total interactions`=0))
          }else{
            d1 <- bind_rows(d1,
                            data.frame(type="Cell-Matrisome",`Percentage of total interactions`=0))
          }
        }
        d1$type <- ifelse(d1$type%in%"Cell-Matrisome","Non.matrisome-Matrisome",d1$type)
        d1$type <- factor(d1$type,levels=c("Non.matrisome-Matrisome","Matrisome-Matrisome"))

        pl <- ggplot(d1,aes(type,`Percentage of total communications`,fill=type)) +
          geom_bar(stat="Identity") +
          #scale_fill_nejm() +
          scale_fill_manual(breaks = d1$type, values = c("#C5F385","#83da07"))+
          ylim(c(0,100))


        pl <- pl + theme_bw()
        # pl <- pl + theme(legend.position = 'none')
        pl <- pl + theme(axis.title.x = element_blank(),
                         axis.text.x = element_blank(),
                         # axis.text.y = element_blank(),
                         axis.ticks.x = element_blank(),
                         panel.grid = element_blank(),
                         plot.title = element_text(size = 20, face = "bold"))
        # pl <- pl + scale_fill_viridis_d(option = "inferno")
        pl <- pl + labs(fill = '') + ggtitle("Communication Pairs")
        return(ggplotly(pl,tooltip = c("Percentage of total communications")))

      }

      #
      #
      #
      #
      #
      #
      #
      # m <- out
      # m <- m[m$`Type of interaction`=="extracellular",]
      # if(nrow(m)<1){
      #   return(ggplot())
      # }else{
      #   m <- distinct(merge(m,mlist,by.x="Gene1",by.y="gene",all.x=T))
      #   m <- distinct(merge(m,mlist,by.x="Gene2",by.y="gene",all.x=T))
      #   r <- as.data.frame(table(c(m$category.x,m$category.y,m$family.x,m$family.y)))
      #   r <- na.omit(r)
      #   dd <- data.frame(Var1=c("Core matrisome",
      #                           "Matrisome-associated",
      #                           "",
      #                           "ECM Glycoproteins",
      #                           "Collagens",
      #                           "Proteoglycans",
      #                           "ECM-affiliated Proteins",
      #                           "ECM Regulators",
      #                           "Secreted Factors"
      #   ),
      #   x=c(1:9))
      #   r <- merge(r,dd,by="Var1",all.y=T)
      #   r[is.na(r)] <- 0
      #   r <- r[order(r$x),]
      #   r$Freq[1] <- round((r$Freq[1]/(r$Freq[1]+r$Freq[2]))*100,1)
      #   r$Freq[2] <- round((r$Freq[2]/(r$Freq[1]+r$Freq[2]))*100,1)
      #   r$Freq[3] <- 0
      #
      #   r$Freq[4] <- round((r$Freq[4]/(r$Freq[4]+r$Freq[5]+r$Freq[6]+r$Freq[7]+r$Freq[8]+r$Freq[9]))*100,1)
      #   r$Freq[5] <- round((r$Freq[5]/(r$Freq[4]+r$Freq[5]+r$Freq[6]+r$Freq[7]+r$Freq[8]+r$Freq[9]))*100,1)
      #   r$Freq[6] <- round((r$Freq[6]/(r$Freq[4]+r$Freq[5]+r$Freq[6]+r$Freq[7]+r$Freq[8]+r$Freq[9]))*100,1)
      #   r$Freq[7] <- round((r$Freq[7]/(r$Freq[4]+r$Freq[5]+r$Freq[6]+r$Freq[7]+r$Freq[8]+r$Freq[9]))*100,1)
      #   r$Freq[8] <- round((r$Freq[8]/(r$Freq[4]+r$Freq[5]+r$Freq[6]+r$Freq[7]+r$Freq[8]+r$Freq[9]))*100,1)
      #   r$Freq[9] <- round((r$Freq[9]/(r$Freq[4]+r$Freq[5]+r$Freq[6]+r$Freq[7]+r$Freq[8]+r$Freq[9]))*100,1)
      #
      #   r <- merge(r, nd1, by.x="Var1",by.y="V1")
      #   colnames(r)[colnames(r)=="Var1"] <- "Matrisome component"
      #   colnames(r)[colnames(r)=="Freq"] <- "Percentage of total interactions"
      #   r$`Matrisome component` <- factor(r$`Matrisome component`,levels=c(
      #     "Core matrisome",
      #     "Matrisome-associated",
      #     "ECM Glycoproteins",
      #     "Collagens",
      #     "Proteoglycans",
      #     "ECM-affiliated Proteins",
      #     "ECM Regulators",
      #     "Secreted Factors"
      #
      #   ))
      #
      #   pl <- ggplot(r,aes(x,`Percentage of total interactions`,fill=`Matrisome component`)) +
      #     geom_bar(stat="Identity")+
      #     scale_fill_manual(breaks = r$`Matrisome component`, values = as.character(r$V2))+
      #     ylim(c(0,100)) +
      #     geom_vline(xintercept = 3, linetype=2, color="grey")
      #
      #   pl <- pl + theme_bw()
      #   # pl <- pl + theme(legend.position = 'none')
      #   pl <- pl + theme(axis.title.x = element_blank(),
      #                    axis.text.x = element_blank(),
      #                    # axis.text.y = element_blank(),
      #                    axis.ticks.x = element_blank(),
      #                    panel.grid = element_blank(),
      #                    plot.title = element_text(size = 20, face = "bold"))
      #  # pl <- pl + scale_fill_viridis_d(option = "inferno")
      #   pl <- pl + labs(fill = '') + ggtitle("Extracellular interactions")
      #   return(ggplotly(pl,tooltip = c("Matrisome component","Percentage of total interactions")))
      # }

    }else{
      return(ggplot())
    }
  }
  output$plot_sub3 <- renderPlotly({
    req(d4())
   # withProgress(message = 'counting the extracellular interactions', value = 0, {
     # print(sub3())
      sub3()
    #  incProgress(1)})
  })

  sub4 <- function(){
    req(d4())
    req(value2()!=0)
    if(attributes(d4())$outcome != "failure"){

      out <- d4()
      out$cc <- paste0(out$Population1," - ",out$Population2)
      out <- out[out$cc %in% selected_state(), ]
      if(nrow(out)<1){
        out <- d4()
        out$cc <- NULL
      }else{
        out$cc <- NULL
      }


      out <- merge(out,mlist,by.x="Gene1",by.y="gene",all.x=T)
      out$family[is.na(out$family)] <- "Non.matrisome"
      out <- merge(out,mlist,by.x="Gene2",by.y="gene",all.x=T)
      out$family.y[is.na(out$family.y)] <- "Non.matrisome"

      df <- as.data.frame(table(out$family.x,out$family.y))
      df <- df[df$Freq>0,]
      if(nrow(df)<1){
        return(ggplot())
      }else{
        tb <- matrix(0,7,7)
        rownames(tb) <- c("ECM glycoproteins",
                          "Collagens",
                          "Proteoglycans",
                          "ECM-affiliated proteins",
                          "ECM regulators",
                          "Secreted factors",
                          "Non.matrisome")

        colnames(tb) <- c("ECM glycoproteins",
                          "Collagens",
                          "Proteoglycans",
                          "ECM-affiliated proteins",
                          "ECM regulators",
                          "Secreted factors",
                          "Non.matrisome")

        for(i in 1:nrow(df)){
          if(tb[rownames(tb)%in%df$Var2[i],
                colnames(tb)%in%df$Var1[i]] == 0
          ){
            tb[rownames(tb)%in%df$Var1[i],
               colnames(tb)%in%df$Var2[i]] <- sum(
                 tb[rownames(tb)%in%df$Var1[i],
                    colnames(tb)%in%df$Var2[i]],
                 df$Freq[i]
               )
          }else{
            tb[rownames(tb)%in%df$Var2[i],
               colnames(tb)%in%df$Var1[i]] <- sum(
                 tb[rownames(tb)%in%df$Var2[i],
                    colnames(tb)%in%df$Var1[i]],
                 df$Freq[i]
               )
          }
        }

        df <- reshape2::melt(tb)
        df <- df[df$value>0,]

        if(nrow(df)<1){
          return(ggplot())
        }else{
          df$variable <- NULL
          ddf <- apply(df,1,function(x){
            if(x[1]%in%"Not matrisome"){
              data.frame(Var1=x[2],Var2=x[1],value=x[3])
            }else{
              data.frame(Var1=x[1],Var2=x[2],value=x[3])
            }
          })
          df <- bind_rows(ddf)
          names(df)[3] <- "Freq"
          df$Freq <- as.numeric(df$Freq)
          df$Freq <- round((df$Freq/sum(df$Freq))*100,3)

          names(df)[1:2] <- c("Gene1","Gene2")

          df$Pair <- paste0(df$Gene1,"-",df$Gene2)
          names(df)[3] <- "Percentage of total communications"
          df <- distinct(merge(df,nd1,by.x="Gene1",by.y="V1",all.x=T))
          df$Gene1 <- factor(df$Gene1,
                             levels = rev(c("ECM glycoproteins",
                                        "Collagens",
                                        "Proteoglycans",
                                        "ECM-affiliated proteins",
                                        "ECM regulators",
                                        "Secreted factors",
                                        "Non.matrisome")))
          df$Gene2 <- factor(df$Gene2,
                             levels = c("ECM glycoproteins",
                                            "Collagens",
                                            "Proteoglycans",
                                            "ECM-affiliated proteins",
                                            "ECM regulators",
                                            "Secreted factors",
                                            "Non.matrisome"))
          pl <- ggplot(df,aes(Gene2,Gene1,label=Pair)) +
            geom_point(aes(size=`Percentage of total communications`,color=as.character(V2))) +
            scale_color_manual(breaks=df$V2,values=as.character(df$V2))
          #+
          # ylim(c(0,100)) +
          # scale_fill_manual(breaks = r$`Matrisome component`, values = as.character(r$V2))+
          # geom_vline(xintercept = 3, linetype=2, color="grey")

          pl <- pl + theme_bw()
          # pl <- pl + theme(legend.position = 'none')
          pl <- pl + theme(# axis.title.x = element_blank(),
                           # axis.title.y = element_blank(),
                           axis.text.x = element_blank(),
                           axis.text.y = element_blank(),
                           axis.ticks.x = element_blank(),
                           axis.ticks.y = element_blank(),
                           panel.grid = element_blank(),
                           plot.title = element_text(size = 20, face = "bold"))
          # pl <- pl + scale_fill_viridis_d(option = "inferno")
          pl <- pl + labs(fill = '') + ggtitle("Matrisome Pairs") + NoLegend()
          return(ggplotly(pl,tooltip = c("Pair","Percentage of total communications")))
          # return(pl)

        }
      }

    #   out <- merge(out,mlist,by.x="Gene1",by.y="gene",all.x=T)
    #   out$family[is.na(out$family)] <- "not matrisome"
    #   out <- merge(out,mlist,by.x="Gene2",by.y="gene",all.x=T)
    #   out$family.y[is.na(out$family.y)] <- "not matrisome"
    #
    #   df <- as.data.frame(table(out$family.x,out$family.y))
    #   df <- df[df$Freq>0,]
    #
    #   if(nrow(df)<1){
    #     return(ggplot())
    #   }else{
    #     colnames(df) <- c("Gene1","Gene2","Freq")
    #     df$Freq <- round((df$Freq/sum(df$Freq))*100,1)
    #     df$Gene1 <- factor(df$Gene1,
    #                        levels = c("ECM Glycoproteins",
    #                                   "Collagens",
    #                                   "Proteoglycans",
    #                                   "ECM-affiliated Proteins",
    #                                   "ECM Regulators",
    #                                   "Secreted Factors",
    #                                   "not matrisome"))
    #     df$Gene2 <- factor(df$Gene2,
    #                        levels = rev(c("ECM Glycoproteins",
    #                                   "Collagens",
    #                                   "Proteoglycans",
    #                                   "ECM-affiliated Proteins",
    #                                   "ECM Regulators",
    #                                   "Secreted Factors",
    #                                   "not matrisome")))
    #     df$pair <- paste0(df$Gene1,"-",df$Gene2)
    #     names(df)[3] <- "percentage of total"
    #     df$pairinv <- paste0(df$Gene2,"-",df$Gene1)
    #     ll <- list()
    #     ex <- list()
    #     for(i in unique(df$pair)){
    #       a <- df[df$pair==i,]
    #       b <- df[df$pair==a$pairinv,]
    #       if(nrow(b)<1){
    #         ll[[i]] <- a
    #       }else{
    #         ex[[i]] <- a$pairinv
    #         if(i %in% unlist(ex)){
    #           next
    #         }else{
    #           z <- a
    #           z$`percentage of total` <- sum(a$`percentage of total`,b$`percentage of total`)
    #           ll[[i]] <- z
    #         }
    #       }
    #     }
    #     ll <- bind_rows(ll)
    #     ll$pair <- NULL
    #     names(ll)[4] <- "pair"
    #     df <- ll
    #     df <- distinct(merge(df,nd1,by.x="Gene2",by.y="V1",all.x=T))
    #
    #     pl <- ggplot(df,aes(Gene1,Gene2,label=pair)) +
    #       geom_point(aes(size=`percentage of total`,color=as.character(V2))) +
    #       scale_color_manual(breaks=df$V2,values=as.character(df$V2))
    #       #+
    #       # ylim(c(0,100)) +
    #       # scale_fill_manual(breaks = r$`Matrisome component`, values = as.character(r$V2))+
    #       # geom_vline(xintercept = 3, linetype=2, color="grey")
    #
    #     pl <- pl + theme_bw()
    #     # pl <- pl + theme(legend.position = 'none')
    #     pl <- pl + theme(axis.title.x = element_blank(),
    #                      axis.title.y = element_blank(),
    #                      axis.text.x = element_blank(),
    #                      axis.text.y = element_blank(),
    #                      axis.ticks.x = element_blank(),
    #                      axis.ticks.y = element_blank(),
    #                      panel.grid = element_blank(),
    #                      plot.title = element_text(size = 20, face = "bold"))
    #     # pl <- pl + scale_fill_viridis_d(option = "inferno")
    #     pl <- pl + labs(fill = '') + ggtitle("Matrisome pairs") + NoLegend()
    #     return(ggplotly(pl,tooltip = c("pair","percentage of total")))
    #   }
    #
    #   # m <- out
    #   # m <- m[m$`Type of interaction`!="extracellular",]
    #   # if(nrow(m)<1){
    #   #   return(ggplot())
    #   # }else{
    #   #   m <- distinct(merge(m,mlist,by.x="Gene1",by.y="gene",all.x=T))
    #   #   m <- distinct(merge(m,mlist,by.x="Gene2",by.y="gene",all.x=T))
    #   #   r <- as.data.frame(table(c(m$category.x,m$category.y,m$family.x,m$family.y)))
    #   #   r <- na.omit(r)
    #   #   dd <- data.frame(Var1=c("Core matrisome",
    #   #                           "Matrisome-associated",
    #   #                           "",
    #   #                           "ECM Glycoproteins",
    #   #                           "Collagens",
    #   #                           "Proteoglycans",
    #   #                           "ECM-affiliated Proteins",
    #   #                           "ECM Regulators",
    #   #                           "Secreted Factors"
    #   #   ),
    #   #   x=c(1:9))
    #   #   r <- merge(r,dd,by="Var1",all.y=T)
    #   #   r[is.na(r)] <- 0
    #   #   r <- r[order(r$x),]
    #   #   r$Freq[1] <- round((r$Freq[1]/(r$Freq[1]+r$Freq[2]))*100,1)
    #   #   r$Freq[2] <- round((r$Freq[2]/(r$Freq[1]+r$Freq[2]))*100,1)
    #   #   r$Freq[3] <- 0
    #   #
    #   #   r$Freq[4] <- round((r$Freq[4]/(r$Freq[4]+r$Freq[5]+r$Freq[6]+r$Freq[7]+r$Freq[8]+r$Freq[9]))*100,1)
    #   #   r$Freq[5] <- round((r$Freq[5]/(r$Freq[4]+r$Freq[5]+r$Freq[6]+r$Freq[7]+r$Freq[8]+r$Freq[9]))*100,1)
    #   #   r$Freq[6] <- round((r$Freq[6]/(r$Freq[4]+r$Freq[5]+r$Freq[6]+r$Freq[7]+r$Freq[8]+r$Freq[9]))*100,1)
    #   #   r$Freq[7] <- round((r$Freq[7]/(r$Freq[4]+r$Freq[5]+r$Freq[6]+r$Freq[7]+r$Freq[8]+r$Freq[9]))*100,1)
    #   #   r$Freq[8] <- round((r$Freq[8]/(r$Freq[4]+r$Freq[5]+r$Freq[6]+r$Freq[7]+r$Freq[8]+r$Freq[9]))*100,1)
    #   #   r$Freq[9] <- round((r$Freq[9]/(r$Freq[4]+r$Freq[5]+r$Freq[6]+r$Freq[7]+r$Freq[8]+r$Freq[9]))*100,1)
    #   #
    #   #   r <- merge(r, nd1, by.x="Var1",by.y="V1")
    #   #   colnames(r)[colnames(r)=="Var1"] <- "Matrisome component"
    #   #   colnames(r)[colnames(r)=="Freq"] <- "Percentage of total interactions"
    #   #
    #   #   r$`Matrisome component` <- factor(r$`Matrisome component`,levels=c(
    #   #     "Core matrisome",
    #   #     "Matrisome-associated",
    #   #     "ECM Glycoproteins",
    #   #     "Collagens",
    #   #     "Proteoglycans",
    #   #     "ECM-affiliated Proteins",
    #   #     "ECM Regulators",
    #   #     "Secreted Factors"
    #   #
    #   #   ))
    #   #
    #   #   pl <- ggplot(r,aes(x,`Percentage of total interactions`,fill=`Matrisome component`)) +
    #   #     geom_bar(stat="Identity")+
    #   #     ylim(c(0,100)) +
    #   #     scale_fill_manual(breaks = r$`Matrisome component`, values = as.character(r$V2))+
    #   #     geom_vline(xintercept = 3, linetype=2, color="grey")
    #   #
    #   #   pl <- pl + theme_bw()
    #   #   # pl <- pl + theme(legend.position = 'none')
    #   #   pl <- pl + theme(axis.title.x = element_blank(),
    #   #                    axis.text.x = element_blank(),
    #   #                    # axis.text.y = element_blank(),
    #   #                    axis.ticks.x = element_blank(),
    #   #                    panel.grid = element_blank(),
    #   #                    plot.title = element_text(size = 20, face = "bold"))
    #   #   # pl <- pl + scale_fill_viridis_d(option = "inferno")
    #   #   pl <- pl + labs(fill = '') + ggtitle("Cell-ECM interactions")
    #   #   return(ggplotly(pl,tooltip = c("Matrisome component","Percentage of total interactions")))
    #   #   }
    #
    }else{
      return(ggplot())
    }
  }
  output$plot_sub4 <- renderPlotly({
    req(d4())
    #withProgress(message = 'counting the cell-ECM interactions', value = 0, {
      sub4()
     # incProgress(1)})
  })

  enplot <- function(){
    req(d4())
    req(value2()!=0)
    if(attributes(d4())$outcome == "failure"){
      return(ggplot())
      }else{

        out <- d4()
        out$cc <- paste0(out$Population1," - ",out$Population2)
        out <- out[out$cc %in% selected_state(), ]
        if(nrow(out)<1){
          out <- d4()
          out$cc <- NULL
        }else{
          out$cc <- NULL
        }

        m <- out
        rr <- matrienrich(m,mlist,signs)
        #rr$alp <- ifelse(rr$p.value<0.05,0.8,0.1)
        s <- unique(rr$signature)
        s1 <- c("Core matrisome","Matrisome-associated",
                "Collagens","Proteoglycans","ECM glycoproteins",
                "ECM regulators","Secreted factors")
        s2 <- s[!(s %in% s1)]
        s2 <- s2[order(s2)]
        s <- c(s1,s2)
        rr$signature <- factor(rr$signature,levels = s)
        #rr$eval <- ifelse(rr$p.value<0.05,NA,rr$p.value)
        if(all(rr$p.value>0.05)){
          return(ggplot())
        }else{
          rr <- rr[rr$p.value<0.05,]
          rr$p.value <- paste0("p value = ",round(rr$p.value,4))
          colnames(rr) <- str_to_title(colnames(rr))

          gp <- ggplot(rr,aes(Signature,Populations,size=Overlap,color=Signature,text=P.value)) +
            geom_point() +
            #scale_y_discrete(position = "right") +
            #geom_point(aes(alpha=alp),show.legend = FALSE) +
            scale_y_discrete(labels = abbreviate) +
            theme_bw() + xlab("") + ylab("") +
            #guides(x =  guide_axis(angle = 90)) +
            theme(axis.title = element_blank(),
                  axis.text.x = element_blank(),
                  #axis.text.y = element_blank(),
                  axis.ticks = element_blank(),
                  panel.grid = element_blank(),
                  plot.title = element_text(size = 20, face = "bold"),
                  plot.title.position = "plot",
                  legend.position = "none")
          gp <- gp + #labs(size = 'Overlap (all p<0.05)') +
            ggtitle("Matrisome-specific Signature Enrichment")
          return(ggplotly(gp,tooltip = c("Populations","color","Overlap","P.value"))%>%
                   layout(
                     title = list(
                       x = 0.01
                     )
                   ))
        }
      }
  }
  output$plot_en <- renderPlotly({
    req(d4())
    #withProgress(message = 'performing matrisome-specific enrichment analysis', value = 0, {
      enplot()
      #incProgress(1)})
  })

  infplot <- function(){
    req(d4())
    req(value2()!=0)
    if(attributes(d4())$outcome == "failure"){
      return(ggplot())
    }else{

      out <- d4()
      out$cc <- paste0(out$Population1," - ",out$Population2)
      out <- out[out$cc %in% selected_state(), ]

      if(nrow(out)<1){
        out <- d4()
        out$cc <- NULL
      }else{
        out$cc <- NULL
      }

      # n <- distinct(data.frame(V1=out$Gene1,V2=out$Gene2))
      n <- distinct(out[,c(2,6)])
      if(nrow(n)<1){
        return(ggplot())
      }else{

        n <- simplify(graph.data.frame(n))

        #va <- page.rank(n)$vector
        va <- degree(n)
        va <- va[order(-va)]
        if(length(va)<=100){va<-va}else{va<-va[1:100]}
        va <- rescale(va,to=c(1,10))
        # vb <- names(neighbors(n,names(va)))
        vb <- lapply(names(va), function(x){
          names(neighbors(n,x))
        })
        vb <- unlist(unique(vb))

        va <- data.frame(V1=names(va),value=va)
        va$x <- c(nrow(va):1)
        l <- list()
        for(i in vb){
          z <- names(neighbors(n,i))
          z <- z[z%in%va$V1]
          if(length(z)<1){next}else{
            nnn <- va$value[va$V1%in%z]/sum(va$value[va$V1%in%z])
            l[[i]] <- data.frame(V1=z,tot=nnn,V2=i)
          }
        }
        l <- bind_rows(l)

        if(nrow(l)<1){
          return(ggplot())
        }else{

        va <- distinct(merge(va,l,by="V1"))
        names(va) <- c("Influencer","size","x","Strength","Influenced")
        va <- distinct(merge(va,mlist,by.x="Influencer",by.y="gene",all.x=T))
        va$category[is.na(va$category)] <- "Non.matrisome"
        va$color <- ifelse(va$category%in%"Core matrisome","#002253",
                           ifelse(va$category%in%"Matrisome-associated","#DB3E18","grey80"))
        names(va)[6] <- "Influencer.Matrisome.Division"
        x <- unique(va$x)[order(unique(va$x))]
        names(x) <- c(1:length(x))
        df <- data.frame(x=x,new.x=names(x))
        va <- distinct(merge(va,df,by="x"))
        va$new.x <- as.factor(va$new.x)

        gp <- ggplot(va,aes(new.x,Influenced,color=color,
                            label=Influencer,
                            # label=Influencer.Matrisome.Division
        )) +
          geom_point(aes(size=Strength)) +
          scale_x_discrete(breaks = unique(va$new.x),labels = unique(va$Influencer)) +
          scale_y_discrete(guide = guide_axis(angle = 90)) +
          scale_color_manual(breaks = as.character(va$color), values = as.character(va$color)) +
          theme_bw() + xlab("Influencers") + ylab("Influenced") +
          theme(#axis.title = element_blank(),
            # axis.text.x = element_blank(),
            # axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank(),
            plot.title = element_text(size = 20, face = "bold"),
            plot.title.position = "plot",
            legend.position = "none") +
          coord_flip()
        gp <- gp + #labs(size = 'Overlap (all p<0.05)') +
          ggtitle("Normalized Influence")


          return(ggplotly(gp,tooltip = c("Influencer","Influenced","Strength"
                                         # ,"Influencer.Matrisome.Division"
                                         ))%>%
                   layout(
                     title = list(
                       x = 0.01
                     )
                   ))
        }
      }




    }
  }
  output$plot_in <- renderPlotly({
    req(d4())
    #withProgress(message = 'performing matrisome-specific enrichment analysis', value = 0, {
    infplot()
    #incProgress(1)})
  })

  observe(
    if(attributes(d2())$outcome != "failure"){

      out <- d4()

      if(nrow(out)<1 | attributes(out)$outcome=="failure"){
        output$tbl2 <- renderDT(data.frame(error.message="No influencers found with these parameters. Please check your selection!"))
      }else{

        out <- d4()
        out$cc <- paste0(out$Population1," - ",out$Population2)
        out <- out[out$cc %in% selected_state(), ]
        if(nrow(out)<1){
          out <- d4()
          out$cc <- NULL
        }else{
          out$cc <- NULL
        }

        n <- distinct(data.frame(V1=out$Gene1,V2=out$Gene2))
        n <- simplify(graph.data.frame(n))

        #va <- page.rank(n)$vector
        va <- degree(n)
        va <- va[order(-va)]
        if(length(va)<=100){va<-va}else{va<-va[1:100]}
        va <- rescale(va,to=c(1,10))
        vb <- names(neighbors(n,names(va)))

        va <- data.frame(V1=names(va),value=va)
        va$x <- c(nrow(va):1)
        l <- list()
        for(i in vb){
          z <- names(neighbors(n,i))
          z <- z[z%in%va$V1]
          if(length(z)<1){next}else{
            nnn <- va$value[va$V1%in%z]/sum(va$value[va$V1%in%z])
            l[[i]] <- data.frame(V1=z,tot=nnn,V2=i)
          }
        }
        l <- bind_rows(l)

        if(length(l)<1){
          output$tbl2 <- renderDT(data.frame(error.message="No interactions found with these parameters. Please check your selection!"))
        }else{
          if(length(va)<1){
            output$tbl2 <- renderDT(data.frame(error.message="No interactions found with these parameters. Please check your selection!"))
          }else{
            va <- distinct(merge(va,l,by="V1"))
            names(va) <- c("Influencer","size","x","strength","Influenced")
            va$strength <- round(va$strength,1)
            names(va)[4] <- "Strength"

            output$tbl2 <- renderDT(DT::datatable(va[,c(1,5,4)], extensions = 'Buttons',
                                                  options = list(#scrollX=TRUE, #lengthMenu = c(5,10,15),
                                                    #paging = TRUE, searching = TRUE,
                                                    #fixedColumns = TRUE, autoWidth = TRUE,
                                                    #ordering = TRUE, dom = 'tB',
                                                    buttons = c('csv', 'excel','pdf')),
                                                  rownames = F)  )
          }
        }




      }
    }else{
      output$tbl2 <- renderDT(data.frame(error.message="No interactions found with these parameters. Please check your selection!")
      )
    }
  ) #table 2

  observe(
    if(attributes(d2())$outcome != "failure"){

      out <- d4()

      if(nrow(out)<1 | attributes(out)$outcome=="failure"){
        output$tbl3 <- renderDT(data.frame(error.message="No enrichments found with these parameters. Please check your selection!"))
      }else{
        out <- d4()
        out$cc <- paste0(out$Population1," - ",out$Population2)
        out <- out[out$cc %in% selected_state(), ]
        if(nrow(out)<1){
          out <- d4()
          out$cc <- NULL
        }else{
          out$cc <- NULL
        }

        m <- out
        rr <- matrienrich(m,mlist,signs)
        #rr$alp <- ifelse(rr$p.value<0.05,0.8,0.1)
        s <- unique(rr$signature)
        s1 <- c("Core matrisome","Matrisome-associated",
                "Collagens","Proteoglycans","ECM glycoproteins",
                "ECM regulators","secreted Factors")
        s2 <- s[!(s %in% s1)]
        s2 <- s2[order(s2)]
        s <- c(s1,s2)
        rr$signature <- factor(rr$signature,levels = s)
        #rr$eval <- ifelse(rr$p.value<0.05,NA,rr$p.value)
        if(all(rr$p.value>0.05)){
          output$tbl3 <- renderDT(data.frame(error.message="No enrichments found with these parameters. Please check your selection!"))
        }else{
          rr <- rr[rr$p.value<0.05,]
          rr$p.value <- paste0("p value = ",round(rr$p.value,4))
          names(rr)[1:3] <- c("Populations","Signature","Overlap")



        output$tbl3 <- renderDT(DT::datatable(rr, extensions = 'Buttons',
                                              options = list(#scrollX=TRUE, #lengthMenu = c(5,10,15),
                                                #paging = TRUE, searching = TRUE,
                                                #fixedColumns = TRUE, autoWidth = TRUE,
                                                #ordering = TRUE, dom = 'tB',
                                                buttons = c('csv', 'excel','pdf')),
                                              rownames = F)  )

        }
      }
    }else{
      output$tbl3 <- renderDT(data.frame(error.message="No enrichments found with these parameters. Please check your selection!")
      )
    }
  ) #table 3


  bbex <- function(){
    req(d4())
    if(attributes(d4())$outcome != "failure"){
      return(savebbplot(d4()))
    }else{
      err <- data.frame(x=1,y=1)
      pl <- ggplot(err,aes(x,y)) + geom_point()
      pl <- pl + theme_bw()
      # pl <- pl + theme(legend.position = 'none')
      pl <- pl + theme(axis.title.x = element_blank(),
                       axis.title.y = element_blank(),
                       axis.text.x = element_blank(),
                       axis.text.y = element_blank(),
                       axis.ticks.x = element_blank(),
                       axis.ticks.y = element_blank(),
                       panel.grid = element_blank(),
                       plot.title = element_text(size = 20, face = "bold"))
      # pl <- pl + scale_fill_viridis_d(option = "inferno")
      pl <- pl + labs(fill = '')
      return(pl)
    }
  }

  sub3ex <- function(){
    req(d4())
    if(attributes(d4())$outcome != "failure"){
      return(savesub3(d4(),selected_state()))
    }else{
      err <- data.frame(x=1,y=1)
      pl <- ggplot(err,aes(x,y)) + geom_point()
      pl <- pl + theme_bw()
      # pl <- pl + theme(legend.position = 'none')
      pl <- pl + theme(axis.title.x = element_blank(),
                       axis.title.y = element_blank(),
                       axis.text.x = element_blank(),
                       axis.text.y = element_blank(),
                       axis.ticks.x = element_blank(),
                       axis.ticks.y = element_blank(),
                       panel.grid = element_blank(),
                       plot.title = element_text(size = 20, face = "bold"))
      # pl <- pl + scale_fill_viridis_d(option = "inferno")
      pl <- pl + labs(fill = '')
      return(pl)
    }
  }

  sub4ex <- function(){
    req(d4())
    if(attributes(d4())$outcome != "failure"){
      return(savesub4(d4(),selected_state(),mlist = mlist))
    }else{
      err <- data.frame(x=1,y=1)
      pl <- ggplot(err,aes(x,y)) + geom_point()
      pl <- pl + theme_bw()
      # pl <- pl + theme(legend.position = 'none')
      pl <- pl + theme(axis.title.x = element_blank(),
                       axis.title.y = element_blank(),
                       axis.text.x = element_blank(),
                       axis.text.y = element_blank(),
                       axis.ticks.x = element_blank(),
                       axis.ticks.y = element_blank(),
                       panel.grid = element_blank(),
                       plot.title = element_text(size = 20, face = "bold"))
      # pl <- pl + scale_fill_viridis_d(option = "inferno")
      pl <- pl + labs(fill = '')
      return(pl)
    }
  }

  plotenex <- function(){
    req(d4())
    if(attributes(d4())$outcome != "failure"){
      return(saveploten(d4(),selected_state(),mlist = mlist,signs = signs))
    }else{
      err <- data.frame(x=1,y=1)
      pl <- ggplot(err,aes(x,y)) + geom_point()
      pl <- pl + theme_bw()
      # pl <- pl + theme(legend.position = 'none')
      pl <- pl + theme(axis.title.x = element_blank(),
                       axis.title.y = element_blank(),
                       axis.text.x = element_blank(),
                       axis.text.y = element_blank(),
                       axis.ticks.x = element_blank(),
                       axis.ticks.y = element_blank(),
                       panel.grid = element_blank(),
                       plot.title = element_text(size = 20, face = "bold"))
      # pl <- pl + scale_fill_viridis_d(option = "inferno")
      pl <- pl + labs(fill = '')
      return(pl)
    }
  }

  plotinex <- function(){
    req(d4())
    if(attributes(d4())$outcome != "failure"){
      return(saveplotin(d4(),selected_state(),mlist))
    }else{
      err <- data.frame(x=1,y=1)
      pl <- ggplot(err,aes(x,y)) + geom_point()
      pl <- pl + theme_bw()
      # pl <- pl + theme(legend.position = 'none')
      pl <- pl + theme(axis.title.x = element_blank(),
                       axis.title.y = element_blank(),
                       axis.text.x = element_blank(),
                       axis.text.y = element_blank(),
                       axis.ticks.x = element_blank(),
                       axis.ticks.y = element_blank(),
                       panel.grid = element_blank(),
                       plot.title = element_text(size = 20, face = "bold"))
      # pl <- pl + scale_fill_viridis_d(option = "inferno")
      pl <- pl + labs(fill = '')
      return(pl)
    }
  }

  # output$mynetworkid <- renderVisNetwork({
  #
  #   withProgress(message = 'creating communication network', value = 0, {
  #
  #   req(d4())
  #   if(attributes(d4())$outcome == "failure"){
  #     return(NULL)
  #   }else{
  #
  #     out <- d4()
  #     out$cc <- paste0(out$Population1," - ",out$Population2)
  #     out <- out[out$cc %in% selected_state(), ]
  #     if(nrow(out)<1){
  #       out <- d4()
  #       out$cc <- NULL
  #     }else{
  #       out$cc <- NULL
  #     }
  #
  #   nodes <- unique(c(out$Gene1,out$Gene2))
  #   nodes <- nodes[order(nodes)]
  #   nodes <- data.frame(id=nodes)
  #   edges <- data.frame(from=out$Gene1,to=out$Gene2)
  #
  #   visNetwork(nodes, edges) %>%
  #     visIgraphLayout() %>%
  #     #visNodes(size = 10) %>%
  #     visOptions(highlightNearest = list(enabled = T, hover = T),
  #                nodesIdSelection = T)  %>%
  #     visInteraction(navigationButtons = TRUE) %>%
  #     visPhysics(stabilization = FALSE) %>%
  #     visEdges(smooth = FALSE)
  #   }
  #     incProgress(1)})
  # })

  # observeEvent(input$img, {
  #   if(input$img=="interaction map"){
  #     ggsave("~/g1.png", sk(), device = "png")
  #     showModal(modalDialog(
  #       title = "interaction map (click anywhere outside the image to close it)",
  #       HTML('<img src="~/g1.png" />'),
  #       easyClose = TRUE,
  #       footer = NULL
  #     ))
  #   }
  # })


  output$export <- downloadHandler( #download all graphs
    filename = function() {
      paste0("MatriCom_graphs_",Sys.time(),".zip")
    },
    content = function(file) {
      showModal(modalDialog("Preparing all files for download and zipping (this might take some time...)", footer=NULL))

      owd <- setwd(tempdir())
      on.exit(setwd(owd))

      png("c1.png",
          width = 4000,
          height = 4000,
          res = 300)
      base::print(bbex())
      dev.off()

      png("c2.png",
          width = 2000,
          height = 2000,
          res = 300)
      base::print(sub3ex())
      dev.off()

      png("c3.png",
          width = 2000,
          height = 2000,
          res = 300)
      base::print(sub4ex())
      dev.off()

      png("c4.png",
          width = 5000,
          height = 5000,
          res = 300)
      base::print(plotinex())
      # print(plotenex())
      dev.off()

      png("c5.png",
          width = 10000,
          height = 10000,
          res = 300)
      # print(plotinex())
      base::print(plotenex())
      dev.off()

      con <- file("README.txt")

      if(is.null(selected_state())){
        writeLines("you have downloaded results for the TOTAL DATASET",con)
      }else{
        writeLines(paste0("you have downloaded results for ", toupper(selected_state())),con)
      }


      zip( file, c("c1.png",
                   "c2.png",
                   "c3.png",
                   "c4.png",
                   "c5.png",
                   "README.txt"
                   ))
      removeModal()
    }
  )

  output$export2 <- downloadHandler( #download the table
    filename = function() {
      paste0("MatriCom_data_",Sys.time(),".XLSX")
    },
    content = function(file) {
      showModal(modalDialog("Preparing tabular data and printing to Excel (this might take some time...)", footer=NULL))
      on.exit(removeModal())

      out <- d4()

      if(nrow(out)<1 | attributes(out)$outcome=="failure"){
        k <- data.frame(error.message="No interactions found with these parameters. Please check your selection!")
        return(k)
      }else{

        out$cc <- paste0(out$Population1," - ",out$Population2)
        out <- out[out$cc %in% selected_state(), ]
        if(nrow(out)<1){
          out <- d4()
          out$cc <- NULL
          out$g1g2 <- NULL
          # out <- out
          out <- distinct(merge(out,mlist,by.x="Gene1",by.y="gene",all.x=T))
          out <- distinct(merge(out,mlist,by.x="Gene2",by.y="gene",all.x=T))
          out[is.na(out)] <- "Non.matrisome"
          names(out)[11:14] <- c("Matrisome.Division.Gene1","Matrisome.Category.Gene1","Matrisome.Division.Gene2","Matrisome.Category.Gene2")
          # out$Reliability.score <- as.numeric(out$Reliability.score)
          # out <- out[order(-out$Reliability.score),]
          # 
          # out <- out[,c(3,2,6,1,9,10,5,8,4,7,11,13)]

        }else{
          out$cc <- NULL
          out$g1g2 <- NULL
          # out <- out
          out <- distinct(merge(out,mlist,by.x="Gene1",by.y="gene",all.x=T))
          out <- distinct(merge(out,mlist,by.x="Gene2",by.y="gene",all.x=T))
          out[is.na(out)] <- "Non.matrisome"
          names(out)[11:14] <- c("Matrisome.Division.Gene1","Matrisome.Category.Gene1","Matrisome.Division.Gene2","Matrisome.Category.Gene2")
          
          # out$Reliability.score <- as.numeric(out$Reliability.score)
          # out <- out[order(-out$Reliability.score),]
          # 
          # out <- out[,c(3,2,6,1,9,10,5,8,4,7,11,13)]
        }

        out$perc.expr.Population1 <- round(out$perc.expr.Population1*100,1)
        out$perc.expr.Population2 <- round(out$perc.expr.Population2*100,1)
        out$mean.expr.Gene1 <- round(out$mean.expr.Gene1,2)
        out$mean.expr.Gene2 <- round(out$mean.expr.Gene2,2)
        out$Reliability.score <- as.numeric(out$Reliability.score)
        out <- out[order(-out$Reliability.score),]
        # out <- out[,c(2:6,1,7:ncol(out))]
        # out <- out[,c(1:11,13)]
        
        out <- out[,c(3,2,6,1,9,10,5,8,4,7,12,11,14,13)]
        out[,3] <- ifelse(out[,3]%in%"cell-matrix","Non.matrisome-Matrisome","Matrisome-Matrisome")
        names(out)[3] <- "Type of communication"
        
        dd <- CCgenes2
        dd$label <- ifelse(dd$label=="extracellular","extracellular (non-matrisome)",
                           ifelse(dd$label=="cell membrane","intracellular",dd$label))
        out <- distinct(merge(out,dd,by.x="Gene1",by.y="gene"))
        out <- distinct(merge(out,dd,by.x="Gene2",by.y="gene"))
        names(out)[15:16] <- c("Compartment.Gene1","Compartment.Gene2")
        out <- out[,c(3,2,4,1,5,6:12,15,13,14,16)]
        
        k <- out

        out <- d4()

        if(nrow(out)<1 | attributes(out)$outcome=="failure"){
          k2 <- data.frame(error.message="No influencers found with these parameters. Please check your selection!")
        }else{

          out <- d4()
          out$cc <- paste0(out$Population1," - ",out$Population2)
          out <- out[out$cc %in% selected_state(), ]
          if(nrow(out)<1){
            out <- d4()
            out$cc <- NULL
          }else{
            out$cc <- NULL
          }

          n <- distinct(data.frame(V1=out$Gene1,V2=out$Gene2))
          n <- simplify(graph.data.frame(n))

          #va <- page.rank(n)$vector
          va <- degree(n)
          va <- va[order(-va)]
          if(length(va)<=100){va<-va}else{va<-va[1:100]}
          va <- rescale(va,to=c(1,10))
          vb <- lapply(names(va), function(x){
            names(neighbors(n,x))
          })
          vb <- unlist(unique(vb))

          va <- data.frame(V1=names(va),value=va)
          va$x <- c(nrow(va):1)
          l <- list()
          for(i in vb){
            z <- names(neighbors(n,i))
            z <- z[z%in%va$V1]
            if(length(z)<1){next}else{
              nnn <- va$value[va$V1%in%z]/sum(va$value[va$V1%in%z])
              l[[i]] <- data.frame(V1=z,tot=nnn,V2=i)
            }
          }
          l <- bind_rows(l)
          if(nrow(l)<1){
            k2 <- data.frame(error.message="No influencers found with these parameters. Please check your selection!")
          }else{
          va <- distinct(merge(va,l,by="V1"))
          names(va) <- c("influencer","size","x","strength","influenced")
          va$strength <- round(va$strength,1)
          names(va)[4] <- "Strength"
          k2 <- va[,c(1,5,4)]
        }}

        out <- d4()

        if(nrow(out)<1 | attributes(out)$outcome=="failure"){
          k3 <- data.frame(error.message="No enrichments found with these parameters. Please check your selection!")
        }else{
          out <- d4()
          out$cc <- paste0(out$Population1," - ",out$Population2)
          out <- out[out$cc %in% selected_state(), ]
          if(nrow(out)<1){
            out <- d4()
            out$cc <- NULL
          }else{
            out$cc <- NULL
          }

          m <- out
          rr <- matrienrich(m,mlist,signs)
          #rr$alp <- ifelse(rr$p.value<0.05,0.8,0.1)
          s <- unique(rr$signature)
          s1 <- c("Core matrisome","Matrisome-associated",
                  "Collagens","Proteoglycans","ECM glycoproteins",
                  "ECM regulators","Secreted factors")
          s2 <- s[!(s %in% s1)]
          s2 <- s2[order(s2)]
          s <- c(s1,s2)
          rr$signature <- factor(rr$signature,levels = s)
          #rr$eval <- ifelse(rr$p.value<0.05,NA,rr$p.value)
          if(all(rr$p.value>0.05)){
            k3 <- data.frame(error.message="No enrichments found with these parameters. Please check your selection!")
          }else{
            rr <- rr[rr$p.value<0.05,]
            rr$p.value <- paste0("p value = ",round(rr$p.value,4))
            names(rr)[1:3] <- c("Populations","Signature","Overlap")
            k3 <- rr

          }
        }

      }


      if(is.null(selected_state())){
        info <- data.frame(status="you have downloaded results for the TOTAL DATASET",
                           file=ifelse(value()!=2,basename(input$file1$name),paste0(input$crt,"_",input$smp)),
                           gene.expression.threshold=input$minexp,
                           positive.population.threshold=input$cellprop,
                           model.maximization=input$mmod,
                           use.exclusion.list=input$excl,
                           reliability.scores=paste(input$postsel,collapse=","),
                           type.of.communication=paste(input$postsel2,collapse=","),
                           cellular.compartments=paste(input$postsel3,collapse=",")
        )
      }else{
        info <- data.frame(status=paste0("you have downloaded results for ", toupper(selected_state())),
                           file=ifelse(value()!=2,basename(input$file1$name),paste0(input$crt,"_",input$smp)),
                           gene.expression.threshold=input$minexp,
                           positive.population.threshold=input$cellprop,
                           model.maximization=input$mmod,
                           use.exclusion.list=input$excl,
                           reliability.scores=paste(input$postsel,collapse=","),
                           type.of.communication=paste(input$postsel2,collapse=","),
                           cellular.compartments=paste(input$postsel3,collapse=",")
        )
      }

      write_xlsx(list("README" = info,"communication network" = k, "network influencers" = k2, "enrichments" = k3), file)


    }
  )


}


# Run the application ----
shinyApp(ui = ui, server = server)
