import numpy as np
import inspect
from .. import teos_air

def attribute_to_function(method_name,function_name,*args):
    #in_args=np.broadcast_arrays(*np.atleast_3d(*args))
    #in_args+=(np.empty_like(in_args[0]),)
    #getattr(teos_air,method_name)(function_name,*in_args)
    #return np.squeeze(in_args[-1])
    return np.squeeze(getattr(teos_air,method_name)(function_name,*np.broadcast_arrays(*np.atleast_3d(*args))))

