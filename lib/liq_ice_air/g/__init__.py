realm=__name__.split('.')[-2]
input_type=__name__.split('.')[-1]

from ... import pyteos_interface
import sys, imp

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

def attfunc(function_name,desc):
    func=lambda A,T,p: pyteos_interface.attribute_to_function(realm+'_'+input_type,function_name,A,T,p)
    func.__name__=function_name
    func.__doc__=docstring.format(desc)
    return func

entropy=                attfunc("entropy","wet air entropy (J/K).")
cond_entropy=           attfunc("cond_entropy","wet air entropy with all the moisture in condensed phase (J/K)")
density=                attfunc("density","wet air density (kg/m3).")
enthalpy=               attfunc("enthalpy","wet air enthalpy (J).")
rh_wmo=                 attfunc("rh_wmo","relative humidity using WMO definition")
temperatureequi=        attfunc("temperatureequi","equivalent temperature (K)")
temperaturewetbulb=     attfunc("temperaturewetbulb","wet bulb temperature (K)")
