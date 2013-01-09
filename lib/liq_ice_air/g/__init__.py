realm=__name__.split('.')[-2]
input_type=__name__.split('.')[-1]

from ... import pyteos_interface
import sys, imp

desc={}
desc['entropy']="wet air entropy (J/K).\n\n\
At the freezing temperature, this function assumes that half of the condensate is in liquid phase and\
that the other half is in ice phase."
desc['cond_entropy']="wet air entropy with all the moisture in condensed phase (J/K)"
desc['rh_wmo']="relative humidity using WMO definition"
desc['temperatureequi']="equivalent temperature (K)"

__all__=desc.keys()

docstring="""
:param A: dry air massfraction (kg/kg)
:type A: np.array.
:param T: absolute temperature (K)
:type T: np.array.
:param p: total pressure (Pa)
:type p: np.array.
:returns: {0}

This is for wet air with ice and liquid.
"""

#Dynamically create function names!
function_module = imp.new_module(realm+'_'+input_type)
function_code="""
def {0}(A,T,p):
    return pyteos_interface.attribute_to_function('{2}','{0}',A,T,p)
{0}.__name__='{0}'
{0}.__doc__=\"\"\"{1}\"\"\"
"""
exec 'from pyteos_air import pyteos_interface' in function_module.__dict__

for function_name in __all__:
    local_function_code=function_code.format(function_name,docstring.format(function_name),realm+'_'+input_type)
    exec local_function_code in function_module.__dict__
    setattr(sys.modules[__name__], function_name, getattr(function_module,function_name))
