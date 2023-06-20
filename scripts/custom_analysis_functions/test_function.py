import pandas as pd
from .AnalysisUtilities import AnalysisUtilities

def custom_analysis(batch_df: pd.DataFrame, uniprot_id: str, PDB_path: str, PAE_path: str):
    
    dataColumns = {'error': []}
    
    AFModel = AnalysisUtilities.load_alphafold(uniprot_id, PDB_path, PAE_path)
    
    if AFModel == False:
        dataColumns['error'] = ['PDB/PAE files not found']*batch_df.shape[0]
    else:
        dataColumns['error'] = ['OK']*batch_df.shape[0]
    
    batch_df = batch_df.assign(**dataColumns)
    
    return (batch_df, uniprot_id)