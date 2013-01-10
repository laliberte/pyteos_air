realm=__name__.split('.')[-2]
input_type=__name__.split('.')[-1]

from ... import pyteos_interface
import sys, imp


docstring="""
:param T: absolute temperature (K)
:type T: np.array.
:param p: total pressure (Pa)
:type p: np.array.
:returns: {0}

This is for wet air with ice and liquid.
"""

def attfunc(function_name,desc):
    func=lambda T,p: pyteos_interface.attribute_to_function(realm+'_'+input_type,function_name,T,p)
    func.__name__=function_name
    func.__doc__=docstring.format(desc)
    return func

massfraction_air=attfunc("massfraction_air","saturation dry air massfraction (kg/kg).")
