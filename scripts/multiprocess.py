from __future__ import annotations
import os
import sys
import multiprocessing as mp
import argparse
import pandas as pd
from time import sleep
from utility_functions import *

argParser = argparse.ArgumentParser(
    prog='analyze_variants.py',
    usage='python multiprocess.py -v <variants filepath> --pdb <pdb dirpath> --pae <json dirpath> [OPTIONS -r/--radius --pae_query_only --plddt_window --scrap_ori --multi]',
    description='Python script to analyze the various aspects of missense variant structural prediction confidence metrics on Alphafold models.'
)

argParser.add_argument('-v', '--variants', help='Absolute/relative filepath to the variants in .csv format.')
argParser.add_argument('-o', '--output', help='Absolute/relative filepath to save the program output.')
argParser.add_argument('--pdb', help='Absolute/relative directory path to the Alphafold structure PDB files.')
argParser.add_argument('--pae', help='Absolute/relative directory path to the Alphafold PAE matrix')
argParser.add_argument('--multi', help='An integer value to specify the number of cores to break up the job into.')
argParser.add_argument('-c', '--custom', help='The name of the module containing custom analysis function. The module MUST EXPORT a function named "custom_analysis". The "custom_analysis" function must take the batched DataFrame and its corresponding uniprot ID as its first two arguments followed by an arbitrary number of positional arguments.')
argParser.add_argument('--args', help='Comma-separated argument string of positional arguments to pass into the custom_analysis function.')

arguments = argParser.parse_args()

# ------------------- GET ARGUMENTS FROM argParser ------------------- #
VARIANTS_PATH = arguments.variants
OUTPUT_PATH = './' if arguments.output == None else arguments.output
PDB_PATH = arguments.pdb
PAE_PATH = arguments.pae
MULTI = 1 if arguments.multi == None else int(arguments.multi)

# custom analysis configuration
CUSTOM_ANALYSIS_MODULE = arguments.custom if arguments.custom != None else None
CUSTOM_ANALYSIS_PACKAGE = 'custom_analysis_functions'

CUSTOM_ANALYSIS_ARGUMENTS = [PDB_PATH, PAE_PATH, *arguments.args.split(',')] if arguments.args != None else [PDB_PATH, PAE_PATH]

taskTracker = mp.Value('i', 0)

def track_progress(result):
    
    with taskTracker.get_lock():
        taskTracker.value -= 1
        print(f'{taskTracker.value} batches left...\r', end='')

# script entrypoint
if __name__ == '__main__':
    
    if MULTI > 1:
        print(f'Multiprocessing enabled. Dispatching work to {MULTI} workers.')
    
    VARIANTS = pd.read_csv(VARIANTS_PATH, sep='\t', header=0)
    
    variantsBatchedByUniprot = batch_by_uniprot(VARIANTS)
    taskQueue = create_task_queue(variantsBatchedByUniprot, *CUSTOM_ANALYSIS_ARGUMENTS)
    
    taskTracker.value = len(taskQueue)
    
    resultsQueue = []
    
    with mp.Pool(initializer=pool_initializer, initargs=('.'.join([CUSTOM_ANALYSIS_PACKAGE, CUSTOM_ANALYSIS_MODULE]),), processes=MULTI) as pool:
        
        for task in taskQueue:
            result = pool.apply_async(worker_func, args=task, callback=track_progress, error_callback=track_progress)
            resultsQueue.append(result)
            
        pool.close()
        pool.join()
        
    resultsQueue = [result.get() for result in resultsQueue]
    resultsQueue.sort(key=lambda x: x[1])
    
    resultDataframes = [result[0] for result in resultsQueue]
    
    
    
    print('program finished')