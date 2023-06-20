import pandas as pd
import importlib

def parse_custom_kwargs(argument_string: str):
    kwargExpressions = argument_string.split(',')
    kwargKeyValueTuples = [kwargExp.split('=') for kwargExp in kwargExpressions]
    kwargKeyValueDict = {kwargTuple[0]: eval(kwargTuple[1]) for kwargTuple in kwargKeyValueTuples}
    return kwargKeyValueDict

# function to break up variants_df into batched by uniprot ID
def batch_by_uniprot(variants_df: pd.DataFrame) -> list:
    '''
    Provide a pandas dataframe to break up into batches by Uniprot ID.
    '''
    
    # group variants into discreet UNIPROT IDs
    variantUniprotIdx: dict[str,list] = {}
    
    for row in variants_df.itertuples():
        uniprot = getattr(row, 'uniprot')
        
        try:
            variantUniprotIdx[uniprot].append(int(getattr(row, 'Index')))
        except KeyError:
            variantUniprotIdx[uniprot] = [int(getattr(row, 'Index'))]
    
    uniprotIds: list[str] = list(variants_df['uniprot'].unique())
    variantsBatchedByUniprot: list[list] = [[variants_df.iloc[variantUniprotIdx[uniprot]], uniprot] for uniprot in uniprotIds]
    
    return variantsBatchedByUniprot


# function to create task queue with positional arguments to pass into worker pool
def create_task_queue(batched_variants: list, *args) -> list:
    '''
    Pass in arguments for the `worker_func` as positional arguments after `batched_variants`. \n
    The returned `taskQueue` has the following signature:\n
    `[(batch1, batch_uniprot1, *arguments for the worker_func), (batch2, batch_uniprot2, *arguments for the worker_func), ...]`
    '''
    taskQueue = [(*batch, *args) for batch in batched_variants]
    return taskQueue


# function to import module into pool processes
def pool_initializer(module_name: str):
    
    # set the custom_analysis function to be global in pool processes
    global custom_analysis
    module = importlib.import_module(module_name)
    custom_analysis = module.custom_analysis


# function to call custom_analysis made global by pool_initializer within each pool process
def worker_func(*task_args):
    return custom_analysis(*task_args)