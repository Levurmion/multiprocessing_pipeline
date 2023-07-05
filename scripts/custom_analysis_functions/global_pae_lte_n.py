import pandas as pd
import numpy as np
from .AnalysisUtilities import AnalysisUtilities as util

def custom_analysis(batch_df: pd.DataFrame, uniprot_id: str, PDB_path: str, PAE_path: str):
    
    dataColumns = {'error': [], 'global_pae_lte_3': [], 'global_pae_lte_4': [], 'global_pae_lte_5': []}
    
    AFModel = util.load_alphafold(uniprot_id, PDB_path, PAE_path)
    
    if AFModel == False:
        dataColumns['error'] = ['PDB/PAE files not found'] * batch_df.shape[0]
        dataColumns['global_pae_lte_3'] = ['PDB/PAE files not found'] * batch_df.shape[0]
        dataColumns['global_pae_lte_4'] = ['PDB/PAE files not found'] * batch_df.shape[0]
        dataColumns['global_pae_lte_5'] = ['PDB/PAE files not found'] * batch_df.shape[0]
    else:
        dataColumns['error'] = ['OK'] * batch_df.shape[0]
        totalPaeEntries = AFModel.get_chain('A').length ** 2
        dataColumns['global_pae_lte_3'] = [np.round(np.sum(AFModel.PAE <= 3)/totalPaeEntries, 5)] * batch_df.shape[0]
        dataColumns['global_pae_lte_4'] = [np.round(np.sum(AFModel.PAE <= 4)/totalPaeEntries, 5)] * batch_df.shape[0]
        dataColumns['global_pae_lte_5'] = [np.round(np.sum(AFModel.PAE <= 5)/totalPaeEntries, 5)] * batch_df.shape[0]
    
    batch_df = batch_df.assign(**dataColumns)
    
    return (batch_df, uniprot_id)