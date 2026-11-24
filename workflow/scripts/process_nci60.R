library(PharmacoGx)
library(dplyr)

rawdata <- Sys.getenv("RAWDATA")
procdata <- Sys.getenv("PROCDATA")
# results <- Sys.getenv("RESULTS")

## TODO: make this run with the config
res.dir <- file.path(procdata, "NCI60")
dir.create(res.dir, showWarnings = FALSE, recursive = TRUE)


pset <- readRDS(file.path(rawdata, "PSets", "NCI60","NCI60_2026.rds"))
pset <- updateObject((pset))


treatment.meta <- treatmentInfo(pset)
treatment.meta <- treatment.meta %>% select("PubChem.CID", "treatmentid", "Compound.Name")
write.csv(treatment.meta,file=file.path(res.dir,"treatment_metadata.csv"))


## STEP 1: We Write The Fit Params A Dataset


response.info <- treatmentResponse(pset)$profiles
sensitivity.info<- treatmentResponse(pset)$sensitivity


print(head(sensitivity.info))
print(head(response.info))
stop()

response.info <- response.info %>%
    select("treatmentid",
    "AnnotationDB.CellLine.Name",
    "Cellosaurus.Accession", "GI50" ,"TGI", "LC50", "IC50",  "HS" , "E0", "E_inf", "EC50","aac_recomputed","ic50_recomputed",
    "Rsq", "ic50_recomputed","aac_recomputed","auc_recomputed","PANEL_CODE", "PANEL_NAME")

keep.tissues <- c("Breast Cancer","Colon Cancer","Small Cell Lung Cancer","Non-Small Cell Lung Cancer")
write.csv(response.info, file = file.path(res.dir,"fit_params.csv"))


for(tissue in keep.tissues){
    temp <- response.info %>% filter(PANEL_NAME==tissue)
    name <- unique(temp$PANEL_CODE)
    
    temp <- temp %>% 
            select("treatmentid","HS" , "E0", "E_inf","IC50",'aac_recomputed') %>%
            group_by(treatmentid) %>%
            summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
            as.data.frame()
    
    temp <- merge(treatment.meta, temp, by="treatmentid") %>% 
            select("treatmentid", "PubChem.CID", "HS" , "E0", "E_inf","IC50",'aac_recomputed')
    temp <- temp[complete.cases(temp),]
    print(head(temp))
    
    write.csv(temp,file.path(res.dir, paste0(name[1],"_params.csv")))


}


sensitivity.info <- sensitivity.info %>% 
        select("treatmentid", "AnnotationDB.CellLine.Name","Cellosaurus.Accession",
        "dose_1", "dose_2","dose_3", "dose_4", "dose_5", "dose_6", "dose_7", "dose_8", "dose_9", "dose_10",
         "viability_1", "viability_2","viability_3", "viability_4", "viability_5" , "viability_6" ,"viability_7", "viability_8", "viability_9", "viability_10",
         "growth_percent_1", "growth_percent_2","growth_percent_3","growth_percent_4", "growth_percent_5", "growth_percent_6", "growth_percent_7", "growth_percent_8","growth_percent_9", "growth_percent_10",
         "PANEL_CODE", "PANEL_NAME","tissueid")

        
write.csv(sensitivity.info, file = file.path(res.dir,"dose_response_curves.csv"))