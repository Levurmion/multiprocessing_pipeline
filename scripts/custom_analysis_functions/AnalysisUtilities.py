from alphafoldmodel import AlphafoldModel
import os

class AnalysisUtilities:
    '''
    A utility class providing static methods for error handling and loading Alphafold models
    '''
    
    def __init__(self):
        pass
    
    @staticmethod
    def check_residues_match(AFmodel: AlphafoldModel, residuePos: int, WTres: str):
        modelResidue = AFmodel.get_residue(residuePos)[0]
        try:
            if modelResidue == WTres:
                return 'match'
            elif modelResidue != WTres:
                return 'no match'
        except ValueError:
            return 'out of range'
    
    @staticmethod
    def append_error_message(error_message: str, data_columns: dict):
        dataColumnNames = list(data_columns.keys())
        for colName in dataColumnNames:
            data_columns[colName].append(error_message)
    
    @staticmethod
    def load_alphafold(uniprotId: str, PDB_path: str, PAE_path: str):
        AF_PDB: str = os.path.join(PDB_path, f'AF-{uniprotId}-F1-model_v4.pdb')
        AF_PAE: str = os.path.join(PAE_path, f'AF-{uniprotId}-F1-predicted_aligned_error_v4.json')
        
        try:
            return AlphafoldModel(AF_PDB, AF_PAE)
        except FileNotFoundError:
            return False