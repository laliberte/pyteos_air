import numpy as np
import inspect
import pyteos_air.teos_air as teos_air

def attribute_to_function(method_name,function_name,*args):
    return np.squeeze(getattr(teos_air,method_name)(function_name,*np.broadcast_arrays(*np.atleast_3d(*args))))

def mp_vec_masked(func,args,pool,num_procs,multi=True):
    #This function simplifies the use of the multiprocessing toolbox.

    dims_frac_proc=np.ma.array(np.array(args[0].shape)/num_procs)
    dim_index=np.argmin(np.ma.masked_where(dims_frac_proc<1.0,dims_frac_proc))
    iter_list=[]
    for in_arr in args:
        iter_list.append(np.split(np.ma.filled(in_arr,fill_value),in_arr.shape[dim_index],axis=dim_index))
    if multi:
        out_var = np.concatenate(pool.map(func,zip(*iter_list)),axis=dim_index)
    else:
        out_var = np.concatenate(map(func,zip(*iter_list)),axis=dim_index)
    return np.ma.masked_where(abs(out_var)>=fill_value,out_var)
