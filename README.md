

# Postprandial Liver - Shiny App
This app displays RNA sequencing data from human liver biopsies under postprandial or fasting conditions for healthy control, NAFLD patients, and cirrhosis patients. 

## Local App Launch
If you wish to view the app online, it is accessible here:
  
https://weweralbrechtsenlab.shinyapps.io/PostprandialLiver/

If you would like to run the app from github, you will need a few R packages:

```{r}
install.packages(c(shiny, EnhancedVolcano, tidyverse, DT, tidyverse, shinyjs, writexl, plotly))
```

Once you have all the packages installed, simply run these lines in R. It will download the app and display it in a browser window:

```{r}
library(EnhancedVolcano)
library(tidyverse)
library(DT)
library(shiny)
library(shinyjs)
library(shinythemes) 
library(writexl)
library(plotly)

runGitHub(rep = "PostprandialLiver_ShinyApp", username = "nicwin98", ref = "main")
```
