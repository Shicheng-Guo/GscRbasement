# Install the rjson package  
#install.packages("rjson")
# Install the rjson package  
#install.packages("foreach")
# Install the devtools package  
#install.packages("devtools") 
# Install rAzureBatch package
#devtools::install_github("Azure/rAzureBatch") 
# Install the doAzureParallel package 
#devtools::install_github("Azure/doAzureParallel") 
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
start_s <- Sys.time() 


end_s <- Sys.time() 

stopCluster(cluster)


