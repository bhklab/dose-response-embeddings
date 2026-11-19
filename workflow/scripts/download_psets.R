library(PharmacoGx)
library(yaml)

rawdata <- Sys.getenv("RAWDATA")
config <- read_yaml(file.path(Sys.getenv("CONFIG"),"download.yaml"))


for(pset in names(config$PSETS))
{
    pset.dir <- file.path(rawdata, "PSets", pset)
    dir.create(pset.dir, showWarnings = FALSE, recursive=TRUE)
    
    print(config$PSETS[[pset]]$name)
    
    pset <- downloadPSet(config$PSETS[[pset]]$name,saveDir = pset.dir)
    pset <- updateObject(pset)
}
