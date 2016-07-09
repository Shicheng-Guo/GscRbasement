install.packages("devtools")
library("devtools")
install_github("ropensci/plotly")
devtools::install_github("ropensci/plotly")

library(plotly)
set_credentials_file("DemoAccount", "lr1c37zw81")
data <- list(
  list(
    x = c(0, 1, 2), 
    y = c(6, 10, 2), 
    error_y = list(
      type = "data", 
      array = c(1, 2, 3), 
      visible = TRUE
    ), 
    type = "scatter"
  )
)
response <- py$plotly(data, kwargs=list(filename="basic-error-bar", fileopt="overwrite"))
url <- response$url
