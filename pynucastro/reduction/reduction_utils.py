
class FailedMPIImport(ImportError):
    
    def __init__(self, *args, **kwargs):
        
        super().__init__(*args, **kwargs)
        
        if self.msg is None:
            self.msg = "Failed to import MPI from mpi4py. Check your mpi4py installation if you" + \
                    " want to run with MPI."
    
    def __getattr__(self, attr):
        
        raise self

def mpi_numpy_decomp(MPI_N, MPI_rank, n):

    if MPI_N <= n[0]:
        
        comp_idx = MPI_rank
        comp_step = MPI_N
        rho_idx = T_idx = 0
        rho_step = T_step = 1
        
    elif MPI_N <= n[0]*n[1]:
        
        m = MPI_N // n[0]
        
        if MPI_rank <= n[0]*m:
            comp_idx = MPI_rank % n[0]
            comp_step = n[0]
            rho_idx = MPI_rank // n[0]
            rho_step = m
        else:
            comp_idx = n[0]
            comp_step = 1
            rho_idx = n[1]
            rho_step = 1
        
        T_idx = 0
        T_step = 1
        
    elif MPI_N <= n[0]*n[1]*n[2]:
        
        m = MPI_N // (n[0] * n[1])
        
        if MPI_rank <= n[0]*n[1]*m:
            comp_idx = MPI_rank % n[0]
            comp_step = n[0]
            rho_idx = (MPI_rank // n[0]) % n[1]
            rho_step = n[1]
            T_idx = MPI_rank // (n[0] * n[1])
            T_step = m
        else:
            comp_idx = n[0]
            comp_step = 1
            rho_idx = n[1]
            rho_step = 1
            T_idx = n[2]
            T_step = 1
            
    else:
        
        m = MPI_N // (n[0] * n[1])
        
        if MPI_rank <= n[0]*n[1]*n[2]:
            comp_idx = MPI_rank % n[0]
            comp_step = n[0]
            rho_idx = (MPI_rank // n[0]) % n[1]
            rho_step = n[1]
            T_idx = MPI_rank // (n[0] * n[1])
            T_step = n[2]
        else:
            comp_idx = n[0]
            comp_step = 1
            rho_idx = n[1]
            rho_step = 1
            T_idx = n[2]
            T_step = 1
            
    return comp_idx, comp_step, rho_idx, rho_step, T_idx, T_step
    
def to_list(x, n=1):
    
    try:
        return list(x)
    except TypeError:
        return [x] * n
