import pandas as pd
import numpy as np
from .AnalysisUtilities import AnalysisUtilities as util

def custom_analysis(batch_df: pd.DataFrame, uniprot_id: str, PDB_path: str, PAE_path: str):
    
    dataColumns = {'uniprot': [], 'position': [], 'error': [], 'plDDT_5.0A': [], 'plDDT_1res': [], 'plDDT_31res': [], 'PAE_5.0A_allPairs': [], 'num_neighbours_5.0A': []}
    
    AFModel = util.load_alphafold(uniprot_id, PDB_path, PAE_path)
    
    if AFModel == False:
        util.append_error_message('PDB/PAE files not found', dataColumns)
        returnDataframe = pd.DataFrame(dataColumns)
        return (returnDataframe, uniprot_id)
    
    modelLength = AFModel.get_chain('A').length
    
    for position in range(1, modelLength+1):
        dataColumns['uniprot'].append(uniprot_id)
        dataColumns['position'].append(position)
        dataColumns['error'].append('OK')
        dataColumns['plDDT_5.0A'].append(AFModel.get_local_plddt(position,radius=5))
        dataColumns['plDDT_1res'].append(AFModel.get_plddt(position)[0])
        dataColumns['plDDT_31res'].append(np.around(AFModel.get_plddt_window(position, window=31)[0],3))
        dataColumns['PAE_5.0A_allPairs'].append(np.around(AFModel.get_local_PAE(position, radius=5, with_query_only=True), 3))
        dataColumns['num_neighbours_5.0A'].append(len(AFModel.get_residues_within(position, radius=5)))
    
    returnDataframe = pd.DataFrame(dataColumns)
    
    return (returnDataframe, uniprot_id)
    
    