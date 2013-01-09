realm=__name__.split('.')[-2]
input_type=__name__.split('.')[-1]

from ... import pyteos_interface
import sys, imp

desc={}
desc['ifl']="isentropic freezing level (Pa).\n\n\
For pressures higher than the **ifl**, there can be only liquid water."
desc['iml']="isentropic melting level (Pa).\n\n\
For pressures lower than the **iml**, there can be only ice."

__all__=desc.keys()

docstring="""
:param A: dry air massfraction (kg/kg)
:type A: np.array.
:param eta: wet air entropy obtained from :func:`pyteos_air.liq_ice_air.g.entropy` (J/K)
:type eta: np.array.
:returns: {0}

This is for wet air with ice and liquid.
"""

#Dynamically create function names!
function_module = imp.new_module(realm+'_'+input_type)
function_code="""
def {0}(A,eta):
    return pyteos_interface.attribute_to_function('{2}','{0}',A,eta)
{0}.__name__='{0}'
{0}.__doc__=\"\"\"{1}\"\"\"
"""
exec 'from pyteos_air import pyteos_interface' in function_module.__dict__

for function_name in __all__:
    local_function_code=function_code.format(function_name,docstring.format(function_name),realm+'_'+input_type)
    exec local_function_code in function_module.__dict__
    setattr(sys.modules[__name__], function_name, getattr(function_module,function_name))
