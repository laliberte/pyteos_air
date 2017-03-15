realm=__name__.split('.')[-2]
input_type=__name__.split('.')[-1]

from ... import pyteos_interface
import sys, imp

docstring="""
:param A: dry air massfraction (kg/kg)
:type A: np.array.
:param eta: wet air entropy obtained from :func:`pyteos_air.liq_ice_air.g.entropy` (J/K)
:type eta: np.array.
:param p: total pressure (Pa)
:type p: np.array.
:returns: {0}

This is for wet air with ice and liquid.
"""

def attfunc(function_name,desc):
    func=lambda A,eta,p: pyteos_interface.attribute_to_function(realm+'_'+input_type,function_name,A,eta,p)
    func.__name__=function_name
    func.__doc__=docstring.format(desc)
    return func

temperature=attfunc("temperature","absolute temperature (K).\n\n\
Note that this function is not an exact inverse of :func:`pyteos_air.liq_ice_air.g.entropy`. \
At the freezing point, this function assumes that a combination of liquid water and ice will \
combine to make sure that the temperature is constant between the Isentropic Freezing level \
:func:`pyteos_air.liq_ice_air.il.ifl` and the Isentropic Melting Level :func:`pyteos_air.liq_ice_air.il.iml`")

density=attfunc("density","density (kg/m3).\n\n\
At the freezing point, this function assumes that a combination of liquid water and ice will \
combine to make sure that the temperature is constant between the Isentropic Freezing level \
:func:`pyteos_air.liq_ice_air.il.ifl` and the Isentropic Melting Level :func:`pyteos_air.liq_ice_air.il.iml`")

enthalpy=attfunc("enthalpy","enthalpy (J).\n\n\
At the freezing point, this function assumes that a combination of liquid water and ice will \
combine to make sure that the temperature is constant between the Isentropic Freezing level \
:func:`pyteos_air.liq_ice_air.il.ifl` and the Isentropic Melting Level :func:`pyteos_air.liq_ice_air.il.iml`")

gc_gv=attfunc("gc_gv","gc_minus_gv (J).\n\n\
At the freezing point, this function assumes that a combination of liquid water and ice will \
combine to make sure that the temperature is constant between the Isentropic Freezing level \
:func:`pyteos_air.liq_ice_air.il.ifl` and the Isentropic Melting Level :func:`pyteos_air.liq_ice_air.il.iml`")

gd_gc=attfunc("gd_gc","gd_minus_gc (J).\n\n\
At the freezing point, this function assumes that a combination of liquid water and ice will \
combine to make sure that the temperature is constant between the Isentropic Freezing level \
:func:`pyteos_air.liq_ice_air.il.ifl` and the Isentropic Melting Level :func:`pyteos_air.liq_ice_air.il.iml`")

gd=attfunc("gd","gd (J).\n\n\
At the freezing point, this function assumes that a combination of liquid water and ice will \
combine to make sure that the temperature is constant between the Isentropic Freezing level \
:func:`pyteos_air.liq_ice_air.il.ifl` and the Isentropic Melting Level :func:`pyteos_air.liq_ice_air.il.iml`")

gc=attfunc("gc","gc (J).\n\n\
At the freezing point, this function assumes that a combination of liquid water and ice will \
combine to make sure that the temperature is constant between the Isentropic Freezing level \
:func:`pyteos_air.liq_ice_air.il.ifl` and the Isentropic Melting Level :func:`pyteos_air.liq_ice_air.il.iml`")

gh_gc=attfunc("gh_gc","gh_minus_gc (J).\n\n\
At the freezing point, this function assumes that a combination of liquid water and ice will \
combine to make sure that the temperature is constant between the Isentropic Freezing level \
:func:`pyteos_air.liq_ice_air.il.ifl` and the Isentropic Melting Level :func:`pyteos_air.liq_ice_air.il.iml`")
