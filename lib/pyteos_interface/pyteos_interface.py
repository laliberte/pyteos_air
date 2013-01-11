import numpy as np
import inspect
import pyteos_air.teos_air as teos_air

def attribute_to_function(method_name,function_name,*args):
    return np.squeeze(getattr(teos_air,method_name)(function_name,*np.broadcast_arrays(*np.atleast_3d(*args))))

