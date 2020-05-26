# young_SHR_rat_artery
## Single-cell RNA sequencing of mesenteric and aortic arteries from young SHR and WKY rats

## 1. Install the latest version of R and Rstudio

https://cran.r-project.org/bin/windows/base/

https://rstudio.com/products/rstudio/download/

## 2. Install packages

Install required packages and dependencies with the following script.

```{r install_packages}
install.packages("shiny")
install.packages("ggplot2")
install.packages("Seurat")
install.packages("cowplot")
install.packages("imager")
install.packages("esquisse")
install.packages("shinyWidgets")
install.packages("DT")
install.packages("plotly")
```

Alternatively, some packages could be installed from github with the following script.

```{r install_packages_github}
install.packages("devtools")
devtools::install_github("rstudio/shiny")
devtools::install_github("tidyverse/ggplot2")
devtools::install_github("dahtah/imager")
devtools::install_github("dreamRs/esquisse")
devtools::install_github("dreamRs/shinyWidgets")
devtools::install_github("rstudio/DT")
devtools::install_github("ropensci/plotly")
```


## 3. Download data from dropbox

Download the folder named "shiny_young_SHR_rat_arteries"

#### https://www.dropbox.com/sh/at7fcag3kpaqme2/AACyjBvBK8v8lRY4RCk3h8SOa?dl=0


## 4. Run the shiny app in Rstudio

```{r run_shinyapp}
shiny::runApp("path/to/shiny_young_SHR_rat_arteries")
# examaple: shiny::runApp("E:/shiny_young_SHR_rat_arteries")

```

## 5. Preview of the web interface

![Image of webpreview](https://raw.githubusercontent.com/WenduoGu/young_SHR_rat_artery/master/Presentation1.jpg)
