# MatriCom (local version), launched by MatriCom::matricom.app().
# Same engine and interface as the online app, without the open-access
# sample collections and without an upload size limit.

suppressPackageStartupMessages({
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
  library(writexl)
  library(packcircles)
  library(ggiraph)
  library(plotly)
  library(viridis)
  library(igraph)
  library(scales)
  library(htmltools)
  library(cicerone)
})

source("matricom_core.R", local = TRUE)
source("io.R", local = TRUE)
source("views.R", local = TRUE)

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

intlist <- readRDS("lr2.RDS")
mlist <- readRDS("matrisome.list.RDS")[[1]]
excl <- as.data.frame(fread("multimers.txt", sep = "\t"))
signs <- as.data.frame(fread("naba-1.csv"))
CCgenes2 <- readRDS("CCgenes2.RDS")
ens2sym <- readRDS("ens2sym.rds") # EnsDb.Hsapiens.v86 GENEID -> SYMBOL

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

matrienrich <- function(data, fams, signs) matrienrich_core(data, fams, signs)

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
  pl <- pl + labs(fill = '') + ggtitle("Matrisome Pairs") + theme(legend.position = "none")
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
  
  n <- simplify(graph_from_data_frame(n))
  
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
                       src = "izzilab logo.png",
                       width = "120px", height = "120px"),
                       href="https://www.oulu.fi/en/research-groups/izzi-group",
                       target="_blank"),tags$a(img(
                         src = "nabalablogo2.png",width = "240px", height = "120px"),
                         href="https://sites.google.com/a/uic.edu/nabalab/",
                         target="_blank")
                     )),
          p(em("If you use MatriCom in your publication, please cite our ",
               a("manuscript",
                 href = "https://doi.org/10.1242/jcs.263927",
                 target="_blank"))),
          
          p("Follow us at ",
            tags$img(
              src = "bluesky.png",
              width = 20,
              height = 20
            ),
            a("The Matrisome Project",
              href = "https://bsky.app/profile/matrisome.bsky.social",
              target="_blank")
          ),
          
          # p("Download the MatrisomeAnalyzeR package ", a("here",
          #                                                href = "https://github.com/Izzilab/MatrisomeAnalyzeR",
          #                                                target="_blank"
          # )),
          
          a(actionButton(inputId = "email1", label = "Contact Us!", 
                         icon = icon("envelope", lib = "font-awesome")),
            href="mailto:matrisomeproject@gmail.com"),
          
          actionButton("guide", "Take a guided tour of MatriCom",
                       icon = icon("circle-info", lib = "font-awesome",
                                   style = "animation: glowing 1300ms infinite;")),
          # Horizontal line
          tags$hr(),
          tags$h2(p("Data Input")),
          tags$h3(p("Upload your dataset")),
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
                                  content = c("MatriCom accepts scRNA-seq files from Seurat or SingleCellExperiment (RDS or QS format), and ScanPy/Loom (H5AD format).",
                                              " ",
                                              "IMPORTANT: MatriCom currently only accepts human and mouse datasets. If your file contains Ensembl gene IDs, you must convert them to HGCN (human) or MGI (mouse) Gene Symbols by activating the conversion button."),
                                  buttonLabel = "OK")),
          
          
          div(style = "margin-top: -17px"),

          div(
            id = "step2",
            awesomeCheckbox(
            inputId = "convbutton",
            label = "Convert my Ensembl IDs to gene symbols", 
            value = FALSE,
            status = "danger"
          )),

          div(style = "margin-top: 20px"),

          tags$hr(),
          
          # Input: Select identities and other params
          tags$h2(p("Query Parameters")),
          div(
            id = "step5",
            uiOutput("id1")
            ),
          div(
            id = "step6",
            sliderInput("minexp", "Set minimum mean gene expression threshold",
                      min = 1,
                      max = 10,
                      step = 1, value = 1) %>% helper(type = "inline",
                                                      title = "Mean gene expression threshold ",
                                                      content = c("This slider defines the minimum mean level above which each gene must be expressed, across the population, for that gene to be included in the analysis."),
                                                      buttonLabel = "OK")
             ),
          div(
            id = "step7",
          sliderInput("cellprop", "Set minimum % positive population threshold",
                        min = 1, max = 100, step = 1, value = 30) %>% helper(type = "inline",
                                                                             title = "Minimum % positive population",
                                                                             content = c("This slider defines the % of cells in a population that must be positive for a gene at the level chosen above for that gene to be included in the analysis."),
                                                                             buttonLabel = "OK")
          ),
          div(
            id = "step10",
            actionButton("runbutton", "Run analysis", icon = icon("upload"), 
                         style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
            
            
            downloadButton("export",
                           label = "Export all graphs (ZIP)",
                           class = "butt"),
            tags$head(tags$style(".butt{background:#cb9d06;} .butt{color: #f5f6f7;}")),
            
            downloadButton("export2", 
                           label = "Export tabular data (XLSX)",
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
                                       content = c("Because the MatriComDB interaction database was curated from multiple sources, analysis may return duplicate entries of the same interaction with differing reliability scores. Selecting the ‘Maximize model’ option will only return interactions with the highest reliability score.",
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
          
          # Horizontal line 
          # tags$hr(),
          
    #       div(
    #         id = "step10",
    #       actionButton("runbutton", "run analysis", icon = icon("upload"), 
    #                    style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
    # 
    # 
    #       downloadButton("export",
    #                      label = "export all graphs (ZIP)",
    #                      class = "butt"),
    #       tags$head(tags$style(".butt{background:#cb9d06;} .butt{color: #f5f6f7;}")),
    #       
    #       downloadButton("export2", 
    #                      label = "export tabular data (XLSX)",
    #                      class = "butt2"),
    #       tags$head(tags$style(".butt2{background:#74b72e;} .butt2{color: #f5f6f7;} .butt2:.noHover{
    # pointer-events: none;}"))),
          
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
            # column(width = 12, withSpinner(echarts4rOutput("plot_gir1",height = 1000
            # ),type=1)),
          ),
          # fluidRow(
          #   column(width = 12, withSpinner(collapsibleTreeOutput("tree",height = 350),type=1)),
          #   
          # ),
          div(style = "height: 50px"),
          fluidRow(
            column(width = 6, withSpinner(plotlyOutput("plot_sub3",height = 350),type=1)),
            column(width = 6, withSpinner(plotlyOutput("plot_sub4",height = 350),type=1)),
          ),
          div(style = "height: 50px"),
          DTOutput('tbl')
            ),
          
          # tabPanel("NETWORK EXPLORER",
          #          # div(style = "height: 20px"),
          #          fluidRow(
          #            # column(width = 12, cyjShinyOutput('cyjShiny'),type=1),
          #            column(width = 12, echarts4rOutput("plot_net",height = 1500),type=1),
          #            # div(style = "height: 50px"),
          #            # DTOutput('tbl3')
          #          ),
          # ),
          
          # tabPanel("NETWORK INFLUENCERS",
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

options(shiny.maxRequestSize = Inf) # no upload limit in the local version

# Define server logic ----
server <- function(input, output, session) {

  observe_helpers()
  guide$init()

  observeEvent(input$guide, {
    guide$start()
  }) #initiates the guide

  selected_state <- reactive({
    input$plot_gir1_selected
  }) #reactive to capture the selection from the bubble plot

  observeEvent(selected_state(), {
    if (!is.null(selected_state())) {
      showNotification(
        paste0("You are now viewing results for ", selected_state(), ".\nYou can remove the selection by clicking it again!"),
        duration = 4, closeButton = T, type = "message")
    }
  }) #side notes for bubble selection

  observeEvent(du(), {
    req(!is.null(du()$data), nrow(du()$data) >= 1)
    rn <- rownames(du()$data)
    nm <- sample(rn, min(10, length(rn)))
    nm <- nm[!(grepl("orf", nm, ignore.case = T))]
    if (all(nm == toupper(nm)) == "FALSE") {
      if (value2() == 56) {value2(1)} else {value2(56)}
      value3(0)
      showModal(modalDialog(
        title = "IMPORTANT!",
        "MatriCom has automatically detected that your input is likely not human. Please note that, in this case, MatriCom directly converts all genes to uppercase before proceeding with the analysis. This should generally work, but if you want to minimize the risk of loosing some genes, please convert your gene IDs before resubmitting your data!",
        easyClose = FALSE,
        footer = modalButton("close and continue")
      ))
    } else {
      if (value2() == 56) {value2(1)} else {value2(56)}
      value3(0)
    }
  }) #updater to value2(); warns about mouse data

  output$id1 <- renderUI({
    selectizeInput("ids", "Select cell identity labels",
                   choices = "",
                   options = list(create = F)) %>% helper(type = "inline",
                                                          title = "select cell identities",
                                                          content = c("Datasets should contain at least one metadata column that holds cell identity labels (e.g., cell type, cluster, etc.) that will be used to define cell populations for analysis.",
                                                                      "Once a file is successfully uploaded, this menu will automatically populate with a list of the column headers."),
                                                          buttonLabel = "OK")
  }) #cell ID button

  #disable all buttons at start - they will pop up when a file is loaded
  for (id in c("minexp","cellprop","scbutton","mmod","excl","nohomo","postsel","postsel2","postsel3","runbutton","export","export2"))
    shinyjs::disable(id)

  value <- reactiveVal(0)   # 0 no data, 1 file loaded
  value2 <- reactiveVal(56) # tracks user/OA switching and governs outputs
  value3 <- reactiveVal(0)

  du <- reactive({
    req(input$file1)
    show_modal_spinner(spin = "self-building-square", color = "grey", text = "Getting your file ready...")
    on.exit(remove_modal_spinner())
    res <- tryCatch(load_input(input$file1$datapath), error = function(e) e)
    if (inherits(res, "error")) {
      showNotification(paste("Could not read the file:", conditionMessage(res)), type = "error", duration = 10)
      value(0)
      req(FALSE)
    }
    if (is.null(res$data) || !("RNA" %in% res$assays)) {
      value(0)
      showNotification("File does not have the expected \"RNA\" slot! You cannot proceed with analysis!", type = "error")
      req(FALSE)
    }
    if (nrow(res$data) < 1) {
      value(0)
      showNotification("File does not have the expected normalized counts! You cannot proceed with analysis!", type = "error")
      req(FALSE)
    }
    value(1)
    showNotification("File correctly read! You can now proceed with analysis!", type = "message")
    res
  }) #user-specific FILE READ-IN

  observeEvent(input$file1, { shinyjs::enable("runbutton") })

  d0 <- reactive({
    if (value() == 1) return(du())
    NULL
  }) #loaded data

  observeEvent(input$runbutton, {
    if (value2() == 56) {value2(1)} else {value2(56)}
  }) #value2 controller

  observe({
    req(!is.null(d0()))
    req(value() != 0 | value2() != 56)
    shinyjs::enable("minexp")
    shinyjs::enable("cellprop")
    shinyjs::enable("scbutton")
  }) #show buttons when the file is correct or changed

  observeEvent(input$file1, { shinyjs::disable("export"); shinyjs::disable("export2") })
  observeEvent(input$runbutton, {
    req(d4())
    shinyjs::enable("export")
    shinyjs::enable("export2")
  }) #activates the download buttons at run

  observe({
    req(value() != 0)
    req(d0())
    nmv <- colnames(d0()$meta)
    updateSelectizeInput(session, "ids", choices = nmv, server = TRUE, options = list(maxOptions = 2000))
  }) #group column name button updater

  observeEvent(input$runbutton, {
    show_modal_spinner(spin = "trinity-rings", color = "blue", text = "Analysis has started...")
    on.exit(remove_modal_spinner())
    req(d4())
    value2(1)
  }) #spinner while the analysis runs

  observe({
    req(d4())
    req(!is_failure(d2()))
    ids <- c("excl","nohomo","postsel","postsel2","postsel3","mmod")
    if (value() != 0 & value2() != 56) for (id in ids) shinyjs::enable(id)
    else for (id in ids) shinyjs::disable(id)
  }) #show the filter buttons when interactions are found

  observe({
    req(d0())
    ids <- c("plot_gir1","plot_in","plot_in2","tree","plot_sub3","plot_sub4","plot_en","tbl","tbl2","tbl3")
    if (value2() == 56) for (id in ids) shinyjs::hide(id)
    else for (id in ids) shinyjs::show(id)
  }) #hide graphs if value2 is zeroed

  observeEvent(input$postsel2, {
    req(d4())
    if (isTruthy(length(input$postsel2) < 1)) {
      showNotification("WARNING: Deselecting both the communication types removes all interactions!", type = "warning")
    }
  })

  observeEvent(input$postsel3, {
    req(d4())
    if (length(intersect("Matrisome", input$postsel3)) < 1) {
      showNotification("WARNING: Deselecting \"Matrisome\" removes all interactions!", type = "error")
    }
  })

  d2 <- eventReactive(input$runbutton, {
    req(d0())
    if (value() == 1) {
      res <- tryCatch({
        df <- matricom_core(d0()$data, d0()$meta[[input$ids]], input$cellprop, input$minexp,
                            intlist, mlist, input$nohomo,
                            conv = isTruthy(input$convbutton), ens2sym = ens2sym)
        if (is_failure(df) || nrow(df) < 1) fail_df() else { attributes(df)$outcome <- "success"; df }
      }, error = function(e) { message("analysis error: ", conditionMessage(e)); fail_df() })
      return(res)
    }
  }) #analysis

  d4 <- reactive({
    req(d2())
    if (!is_failure(d2())) {
      filtres_core(d2(), input$postsel, input$excl, input$postsel2, input$nohomo, input$postsel3,
                   dups = input$mmod, excl = excl, CCgenes2 = CCgenes2)
    } else d2()
  }) #post-run filters

  # rows for the current bubble selection (all rows if nothing selected)
  sel_out <- reactive({
    req(d4())
    out <- d4()
    if (is_failure(out)) return(out)
    cc <- paste0(out$Population1, " - ", out$Population2)
    keep <- cc %in% selected_state()
    if (any(keep)) out[keep, ] else out
  })

  enr <- reactive({
    out <- sel_out()
    req(!is_failure(out))
    rr <- matrienrich(out, mlist, signs)
    rr
  })

  infl <- reactive({
    out <- sel_out()
    req(!is_failure(out))
    influencers_core(out)
  })

  output$tbl <- renderDT({
    req(d2())
    if (is_failure(d2())) return(datatable(fail_df()))
    out <- d4()
    if (nrow(out) < 1 || is_failure(out)) return(datatable(fail_df()))
    out <- comm_table(sel_out(), mlist)
    DT::datatable(out, extensions = 'Buttons',
                  options = list(buttons = c('csv', 'excel', 'pdf')), rownames = F) %>%
      formatStyle('Reliability.score', background = styleColorBar(c(0, 3), 'lightblue'),
                  backgroundSize = '95% 80%', backgroundRepeat = 'no-repeat', backgroundPosition = 'left') %>%
      formatStyle('perc.expr.Population1', background = styleColorBar(c(0, 100), 'lightcoral'),
                  backgroundSize = '95% 80%', backgroundRepeat = 'no-repeat', backgroundPosition = 'left') %>%
      formatStyle('perc.expr.Population2', background = styleColorBar(c(0, 100), 'lightcoral'),
                  backgroundSize = '95% 80%', backgroundRepeat = 'no-repeat', backgroundPosition = 'left') %>%
      formatStyle('mean.expr.Gene1', background = styleColorBar(c(0, max(out$mean.expr.Gene1)), 'bisque'),
                  backgroundSize = '95% 80%', backgroundRepeat = 'no-repeat', backgroundPosition = 'left') %>%
      formatStyle('mean.expr.Gene2', background = styleColorBar(c(0, max(out$mean.expr.Gene2)), 'bisque'),
                  backgroundSize = '95% 80%', backgroundRepeat = 'no-repeat', backgroundPosition = 'left')
  }) #table

  output$plot_gir1 <- renderGirafe({
    req(d4())
    req(value2() != 0)
    if (!is_failure(d4()) && nrow(d4()) >= 1) bbplot(d4()) else NULL
  }) #bubble plot

  output$plot_sub3 <- renderPlotly({
    out <- sel_out()
    req(value2() != 0)
    if (is_failure(out)) return(ggplotly(ggplot()))
    sub3_plot(out)
  })

  output$plot_sub4 <- renderPlotly({
    out <- sel_out()
    req(value2() != 0)
    if (is_failure(out)) return(ggplotly(ggplot()))
    sub4_plot(out, mlist)
  })

  output$plot_en <- renderPlotly({
    out <- sel_out()
    req(value2() != 0)
    if (is_failure(out)) return(ggplotly(ggplot()))
    enrich_plot(enr())
  })

  output$plot_in <- renderPlotly({
    out <- sel_out()
    req(value2() != 0)
    if (is_failure(out)) return(ggplotly(ggplot()))
    infl_plot(infl(), mlist)
  })

  output$tbl2 <- renderDT({
    req(d2())
    msg <- "No interactions found with these parameters. Please check your selection!"
    if (is_failure(d2())) return(datatable(data.frame(error.message = msg)))
    out <- d4()
    if (nrow(out) < 1 || is_failure(out))
      return(datatable(data.frame(error.message = "No influencers found with these parameters. Please check your selection!")))
    va <- infl()
    if (is.null(va)) return(datatable(data.frame(error.message = msg)))
    names(va) <- c("Influencer", "size", "x", "strength", "Influenced")
    va$strength <- round(va$strength, 1)
    names(va)[4] <- "Strength"
    DT::datatable(va[, c(1, 5, 4)], extensions = 'Buttons',
                  options = list(buttons = c('csv', 'excel', 'pdf')), rownames = F)
  }) #table 2

  output$tbl3 <- renderDT({
    req(d2())
    msg <- "No enrichments found with these parameters. Please check your selection!"
    if (is_failure(d2())) return(datatable(data.frame(error.message = msg)))
    out <- d4()
    if (nrow(out) < 1 || is_failure(out)) return(datatable(data.frame(error.message = msg)))
    rr <- enrich_table(enr())
    if (is.null(rr)) return(datatable(data.frame(error.message = msg)))
    DT::datatable(rr, extensions = 'Buttons', options = list(buttons = c('csv', 'excel', 'pdf')), rownames = F)
  }) #table 3

  output$export <- downloadHandler( #download all graphs
    filename = function() {
      paste0("MatriCom_graphs_", Sys.time(), ".zip")
    },
    content = function(file) {
      showModal(modalDialog("Preparing all files for download and zipping (this might take some time...)", footer = NULL))
      on.exit(removeModal(), add = TRUE)
      d <- tempfile("mcexp"); dir.create(d)
      ok <- !is_failure(d4())
      pngs <- list(
        c1.png = list(if (ok) savebbplot(d4()) else blank_plot(), 4000, 4000),
        c2.png = list(if (ok) savesub3(d4(), selected_state()) else blank_plot(), 2000, 2000),
        c3.png = list(if (ok) savesub4(d4(), selected_state(), mlist = mlist) else blank_plot(), 2000, 2000),
        c4.png = list(if (ok) saveplotin(d4(), selected_state(), mlist) else blank_plot(), 5000, 5000),
        c5.png = list(if (ok) saveploten(d4(), selected_state(), mlist = mlist, signs = signs) else blank_plot(), 10000, 10000)
      )
      for (f in names(pngs)) {
        png_device(file.path(d, f), width = pngs[[f]][[2]], height = pngs[[f]][[3]], res = 300)
        base::print(pngs[[f]][[1]])
        dev.off()
      }
      writeLines(if (is.null(selected_state())) "you have downloaded results for the TOTAL DATASET"
                 else paste0("you have downloaded results for ", toupper(selected_state())),
                 file.path(d, "README.txt"))
      zip::zip(file, files = c(names(pngs), "README.txt"), root = d)
    }
  )

  output$export2 <- downloadHandler( #download the table
    filename = function() {
      paste0("MatriCom_data_", Sys.time(), ".XLSX")
    },
    content = function(file) {
      showModal(modalDialog("Preparing tabular data and printing to Excel (this might take some time...)", footer = NULL))
      on.exit(removeModal())
      out <- d4()
      if (nrow(out) < 1 || is_failure(out)) {
        k <- k2 <- k3 <- data.frame(error.message = "No interactions found with these parameters. Please check your selection!")
      } else {
        k <- comm_export(sel_out(), mlist, CCgenes2)
        va <- infl()
        if (is.null(va)) {
          k2 <- data.frame(error.message = "No influencers found with these parameters. Please check your selection!")
        } else {
          names(va) <- c("influencer", "size", "x", "strength", "influenced")
          va$strength <- round(va$strength, 1)
          names(va)[4] <- "Strength"
          k2 <- va[, c(1, 5, 4)]
        }
        k3 <- enrich_table(enr())
        if (is.null(k3)) k3 <- data.frame(error.message = "No enrichments found with these parameters. Please check your selection!")
      }
      info <- data.frame(status = if (is.null(selected_state())) "you have downloaded results for the TOTAL DATASET"
                                  else paste0("you have downloaded results for ", toupper(selected_state())),
                         file = basename(input$file1$name),
                         gene.expression.threshold = input$minexp,
                         positive.population.threshold = input$cellprop,
                         model.maximization = input$mmod,
                         use.exclusion.list = input$excl,
                         reliability.scores = paste(input$postsel, collapse = ","),
                         type.of.communication = paste(input$postsel2, collapse = ","),
                         cellular.compartments = paste(input$postsel3, collapse = ","))
      write_xlsx(list("README" = info, "communication network" = k, "network influencers" = k2, "enrichments" = k3), file)
    }
  )
}

# Run the application ----
shinyApp(ui = ui, server = server)
