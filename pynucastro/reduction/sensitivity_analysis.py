import numpy as np
from mpi4py import MPI

def binary_search(network, nuclei, errfunc, thresh=0.05):
    
    start_idx = 0
    seg_size = len(nuclei) / 2
    
    while seg_size >= 0.5:
            
        # Divide up into segments
        end_idx = start_idx + round(seg_size)
        red_net = network.linking_nuclei(nuclei[:end_idx])
                
        # Evaluate error
        err = errfunc(red_net)
        
        if err <= thresh:
            seg_size /= 2
        else:
            start_idx += round(seg_size)
            seg_size /= 2
    
    return network.linking_nuclei(nuclei[:start_idx+1])

def n_ary_search(network, nuclei, errfunc, thresh=0.05):
    
    comm = MPI.COMM_WORLD
    MPI_N = comm.Get_size()
    MPI_rank = comm.Get_rank()
    
    start_idx = 0
    seg_size = len(nuclei) / (MPI_N + 1)
    
    while True:
        
        procs_needed = round((MPI_N + 1) * seg_size)
        
        if MPI_rank < procs_needed:
            
            # Divide up into segments
            end_idx = start_idx + round((MPI_rank + 1) * seg_size)
            red_net = network.linking_nuclei(nuclei[:end_idx])
                
            # Evaluate error
            err = errfunc(red_net)
            locflag = np.array([err <= thresh], dtype=np.bool_)
            
        else:
            
            locflag = np.array([True], dtype=np.bool_)

        # Transmit verdict on error
        flags = np.zeros(MPI_N, dtype=np.bool_)
        comm.Allgather([locflag, MPI.BOOL], [flags, MPI.BOOL])
        
        # Figure out what region to search next
        i1 = flags.argmax()
        
        if flags[i1]:
            start_idx += round(seg_size * i1)
            seg_size  /= (MPI_N + 1)
        else:
            start_idx += round(MPI_N * seg_size)
            seg_size  /= (MPI_N + 1)
        
        if procs_needed <= MPI_N:
            break
    
    return network.linking_nuclei(nuclei[:start_idx+1])
