# Install the rjson package  
 install.packages("rjson")

# Install the rjson package  
 install.packages("foreach")

# Install the devtools package  
 install.packages("devtools") 

# Install rAzureBatch package
 devtools::install_github("Azure/rAzureBatch") 

# Install the doAzureParallel package 
devtools::install_github("Azure/doAzureParallel") 

# Load the doAzureParallel library 
library(doAzureParallel) 
#generateCredentialsConfig("credentials.json") 
#generateClusterConfig("cluster.json")
setCredentials("credentials.json")

# Create your cluster if it does not exist; this takes a few minutes
cluster <- makeCluster("cluster.json") 

# Register your parallel backend 
registerDoAzureParallel(cluster) 

# Check that the nodes are running 
getDoParWorkers() 

mean_change = 1.001 
volatility = 0.01 
opening_price = 100 

getClosingPrice <- function() { 
  days <- 1825 # ~ 5 years 
  movement <- rnorm(days, mean=mean_change, sd=volatility) 
  path <- cumprod(c(opening_price, movement)) 
  closingPrice <- path[days] 
  return(closingPrice) 
} 

start_s <- Sys.time() 
# Run 10,000 simulations in series 
closingPrices_s <- foreach(i = 1:10, .combine='c') %do% { 
  replicate(10000, getClosingPrice()) 
} 
write.table(closingPrices_s,file="x.txt")
end_s <- Sys.time() 

hist(closingPrices_s)
stopCluster(cluster)


