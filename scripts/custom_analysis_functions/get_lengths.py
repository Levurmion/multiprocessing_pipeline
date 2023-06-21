import pandas as pd
from .AnalysisUtilities import AnalysisUtilities as util

def custom_analysis(batch_df: pd.DataFrame, uniprot_id: str, PDB_path: str, PAE_path: str):
    
    dataColumns = {'error': [], 'length': []}
    
    AFModel = util.load_alphafold(uniprot_id, PDB_path, PAE_path)
    
    if AFModel == False:
        dataColumns['error'] = ['PDB/PAE files not found']
        dataColumns['length'] = ['PDB/PAE files not found']
    else:
        dataColumns['error'] = ['OK']
        dataColumns['length'] = [AFModel.get_chain('A').length]
    
    batch_df = batch_df.assign(**dataColumns)
    
    return (batch_df, uniprot_id)