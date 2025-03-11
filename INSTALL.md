## Local MatriCom Installation
The online MatriCom app has an upload limit of 1 GB, which may not be suitable for more sizeable files. In this case, we suggest installing MatriCom locally. Below is the list of all required packages (MatriCom is last) and how to install them one by one, if necessary:

```r
# from CRAN:
install.packages("cicerone")
install.packages("shiny")
install.packages("Matrix")
install.packages("DT")
install.packages("dplyr")
install.packages("data.table")
install.packages("ggplot2")
install.packages("shinyWidgets")
install.packages("shinybusy")
install.packages("shinyjs")
install.packages("shinyhelper")
install.packages("shinycssloaders")
install.packages("shinyBS")
install.packages("Seurat")
install.packages("SeuratObject")
install.packages("scCustomize")
install.packages("sp")
install.packages("uwot")
install.packages("shinydashboard")
install.packages("gridExtra")
install.packages("writexl")
install.packages("packcircles")
install.packages("ggiraph")
install.packages("plotly")
install.packages("viridis")
install.packages("collapsibleTree")
install.packages("igraph")
install.packages("scales")
install.packages("scater")
install.packages("qs")
install.packages("stringr")
install.packages("reshape2")
install.packages("htmltools")

# from BioConductor:
install.packages("BiocManager")
BiocManager::install("ensembldb")
BiocManager::install("EnsDb.Hsapiens.v86")

# from GitHub:
devtools::install_github("mojaveazure/seurat-disk")
devtools::install_github("izzilab/MatriCom")

# when ready, run:
matricom.app()
```
