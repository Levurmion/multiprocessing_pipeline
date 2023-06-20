from .AnalysisUtilities import AnalysisUtilities as util
import pandas as pd
import numpy as np

def custom_analysis(batch_df: pd.DataFrame, uniprot_id: str, PDB_path: str, PAE_path: str, window_size: int):
    
    window_size_int = int(window_size)
    windowColname = f'highest_mean_plDDT_in_{window_size}res'
    
    dataColumns = {'error': [], windowColname: []}
    
    AFModel = util.load_alphafold(uniprot_id, PDB_path, PAE_path)
    
    if AFModel == False:
        dataColumns['error'] = ['PDB/PAE files not found']*batch_df.shape[0]
        dataColumns[windowColname] = ['PDB/PAE files not found']*batch_df.shape[0]
        batch_df = batch_df.assign(**dataColumns)
        
        return (batch_df, uniprot_id)
    
    for row in batch_df.itertuples():
        position = getattr(row, 'position')
        WTaa = getattr(row, 'WT')
        
        WTMatches = util.check_residues_match(AFModel, position, WTaa)
        
        if WTMatches == 'match':
            
            positionShift = window_size_int//2
            scanStartPosition = position - positionShift if position - positionShift > 1 else 1
            scanEndPosition = position + positionShift if position + positionShift < AFModel.get_chain('A').length else AFModel.get_chain('A').length
            scanIntervalPlddt = []
            
            for scanPosition in range(scanStartPosition, scanEndPosition + 1):
                intervalPlddt = AFModel.get_plddt_window(scanPosition, window=window_size_int)[0]
                scanIntervalPlddt.append(np.around(intervalPlddt, 3))
            
            dataColumns['error'].append('OK')
            dataColumns[windowColname].append(max(scanIntervalPlddt))
            
        elif WTMatches == 'no match':
            util.append_error_message('WT residue does not match', dataColumns)
            
        elif WTMatches == 'out of range':
            util.append_error_message('query residue out of range', dataColumns)

    batch_df = batch_df.assign(**dataColumns)
    
    return (batch_df, uniprot_id)
