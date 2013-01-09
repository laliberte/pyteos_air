realm=__name__.split('.')[-2]
input_type=__name__.split('.')[-1]

from ... import pyteos_interface
import sys, imp

desc={}
desc['pottemp']="potential tempeature (K) at **pref**."
desc['pottempequisat']="saturation equivalent potential temperature (K) at **pref**."
desc['pottempequipseudo']="pseudo equivalent potential temperature (K) at **pref**.\n\n\
Integrates the pseudo adiabat until the temperature falls\
outside the range of validity of TEOS-10."
desc['pottempequi']="equivalent potential tempeature (K) at **pref**."

__all__=desc.keys()

docstring="""
:param A: dry air massfraction (kg/kg)
:type A: np.array.
:param T: absolute temperature (K)
:type T: np.array.
:param p: total pressure (Pa)
:type p: np.array.
:param pref: reference total pressure (Pa)
:type pref: np.array.
:returns: {0}

This is for wet air with ice and liquid.
"""

#Dynamically create function names!
function_module = imp.new_module(realm+'_'+input_type)
function_code="""
def {0}(A,T,p,pref):
    return pyteos_interface.attribute_to_function('{2}','{0}',A,T,p,pref)
{0}.__name__='{0}'
{0}.__doc__=\"\"\"{1}\"\"\"
"""
exec 'from pyteos_air import pyteos_interface' in function_module.__dict__

for function_name in __all__:
    local_function_code=function_code.format(function_name,docstring.format(function_name),realm+'_'+input_type)
    exec local_function_code in function_module.__dict__
    setattr(sys.modules[__name__], function_name, getattr(function_module,function_name))
