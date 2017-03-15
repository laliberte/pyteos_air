realm=__name__.split('.')[-2]
input_type=__name__.split('.')[-1]

from ... import pyteos_interface
import sys, imp

docstring="""
:param A: dry air massfraction (kg/kg)
:type A: np.array.
:param eta: wet air entropy obtained from :func:`pyteos_air.liq_ice_air.g.entropy` (J/K)
:type eta: np.array.
:returns: {0}

This is for wet air with ice and liquid.
"""

def attfunc(function_name,desc):
    func=lambda A,eta: pyteos_interface.attribute_to_function(realm+'_'+input_type,function_name,A,eta)
    func.__name__=function_name
    func.__doc__=docstring.format(desc)
    return func

ifl=        attfunc("ifl","isentropic freezing level (Pa).\n\n\
For pressures higher than the **ifl**, there can be only liquid water.")
iml=        attfunc("iml","isentropic melting level (Pa).\n\n\
For pressures lower than the **iml**, there can be only ice.")
