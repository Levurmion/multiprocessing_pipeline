import pandas as pd
from .AnalysisUtilities import AnalysisUtilities

def custom_analysis(batch_df: pd.DataFrame, uniprot_id: str, PDB_path: str, PAE_path: str):
    
    AFModel = AnalysisUtilities.load_alphafold(uniprot_id, PDB_path, PAE_path)
    
    if AFModel == False:
        print('no model!')
    else:
        print('model found!')
    
    return (batch_df, uniprot_id)