import sys
import numpy as np 
import pandas as pd
import iisignature as iis
from damply import dirs
from collections import defaultdict
import tqdm


data = pd.read_csv(dirs.PROCDATA / "NCI60/dose_response_curves.csv",index_col=0)



BASE_POINT = np.array([np.amin(data['dose_1']),1])
MAX_LEVELS = 5

treatment_info = pd.read_csv(dirs.PROCDATA /"NCI60/treatment_metadata.csv",index_col=0)
treatment_info = treatment_info.dropna(subset=['PubChem.CID'])
treatment_info = treatment_info.drop_duplicates()


treatment_map = {
    treatment_info['treatmentid'].values[j]: int(treatment_info['PubChem.CID'].values[j]) 
    for j in range(treatment_info.shape[0])
}

signatures = defaultdict(list)
logsignatures = defaultdict(list)

growth_signatures = defaultdict(list)
growth_logsignatures = defaultdict(list)

for drug in (pbar:=tqdm.tqdm(pd.unique(data['treatmentid']))):
    pbar.set_description(f"Working on: {drug}")

    if drug not in treatment_map:
        continue

    temp =data[data['treatmentid']==drug]
    
    for idx, row in tqdm.tqdm(temp.iterrows(),total=temp.shape[0],leave=False):
        # print(row)
        # sys.exit()
        doses = [row[x] for x in row.index if x.startswith('dose')]
        viabilities = [row[x] for x in row.index if x.startswith('viability')]
        growth =  [row[x] for x in row.index if x.startswith('growth')]
        # print(doses)
        doses = [x for x in doses if not np.isnan(x)]
        viabilities = [x/100 for x in viabilities if not np.isnan(x)]
        growth = [x/100 for x in growth if not np.isnan(x)]

        viability_mismatch = len(doses)!=len(viabilities)
        growth_mismatch = len(doses)!=len(growth)

        if not viability_mismatch:
            
            viability_curve = np.vstack( (
                    np.array(doses),
                    np.array(viabilities)
            )).T

            viability_curve = np.vstack((BASE_POINT,viability_curve))
            
            sig = iis.sig(viability_curve ,3)
            s = iis.prepare(2 ,3,"S")
            logsig = iis.logsig(viability_curve,s)

            signatures['treatmentid'].append(drug)
            signatures['PubChem.CID'].append(treatment_map[row['treatmentid']])
            signatures['Cell Line Name'].append(row["AnnotationDB.CellLine.Name"])
            signatures['PANEL_CODE'].append(row['PANEL_CODE'])
            signatures['PANEL_NAME'].append(row['PANEL_NAME'])
            signatures['Tissue'].append(row['tissueid'])
            for j in range(len(sig)):
                signatures[f"Emb-{j+1}"].append(sig[j])



            logsignatures['treatmentid'].append(drug)
            logsignatures['PubChem.CID'].append(treatment_map[row['treatmentid']])
            logsignatures['Cell Line Name'].append(row["AnnotationDB.CellLine.Name"])
            logsignatures['PANEL_CODE'].append(row['PANEL_CODE'])
            logsignatures['PANEL_NAME'].append(row['PANEL_NAME'])
            logsignatures['Tissue'].append(row['tissueid'])
            for j in range(len(logsig)):
                logsignatures[f"Emb-{j+1}"].append(logsig[j])

    # if not growth_mismatch:
                
    #     growth_curve = np.vstack( (
    #             np.array(doses),
    #             np.array(growth)
    #     )).T
        

    #     # viability_curve = np.vstack((BASE_POINT,viability_curve))
        
    #     sig = iis.sig(growth_curve ,3)
    #     s = iis.prepare(2 ,3,"S")
    #     logsig = iis.logsig(growth_curve,s)

    #     growth_signatures['treatmentid'].append(drug)
    #     growth_signatures['PubChem.CID'].append(treatment_map[row['treatmentid']])
    #     growth_signatures['Cell Line Name'].append(row["AnnotationDB.CellLine.Name"])
    #     growth_signatures['PANEL_CODE'].append(row['PANEL_CODE'])
    #     growth_signatures['PANEL_NAME'].append(row['PANEL_NAME'])
    #     growth_signatures['Tissue'].append(row['tissueid'])
    #     for j in range(len(sig)):
    #         growth_signatures[f"Emb-{j+1}"].append(sig[j])



    #     growth_logsignatures['treatmentid'].append(drug)
    #     growth_logsignatures['PubChem.CID'].append(treatment_map[row['treatmentid']])
    #     growth_logsignatures['Cell Line Name'].append(row["AnnotationDB.CellLine.Name"])
    #     growth_logsignatures['PANEL_CODE'].append(row['PANEL_CODE'])
    #     growth_logsignatures['PANEL_NAME'].append(row['PANEL_NAME'])
    #     growth_logsignatures['Tissue'].append(row['tissueid'])
    #     for j in range(len(logsig)):
    #         growth_logsignatures[f"Emb-{j+1}"].append(logsig[j])

signatures = pd.DataFrame(signatures)
logsignatures =pd.DataFrame(logsignatures)

# growth_signatures = pd.DataFrame(growth_signatures)
# growth_logsignatures = pd.DataFrame(growth_logsignatures)

signatures.to_csv(dirs.RESULTS / "viability_signatures.csv")
logsignatures.to_csv(dirs.RESULTS / "viability_log_signatures.csv")

# growth_signatures.to_csv(dirs.RESULTS / "growth_signatures.csv")
# growth_logsignatures.to_csv(dirs.RESULTS / "growth_log_signatures.csv")


