realm=__name__.split('.')[-2]
input_type=__name__.split('.')[-1]

from ... import pyteos_interface

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

def attfunc(function_name,desc):
    func=lambda A,T,p,pref: pyteos_interface.attribute_to_function(realm+'_'+input_type,function_name,A,T,p,pref)
    func.__name__=function_name
    func.__doc__=docstring.format(desc)
    return func

#Functions list:
pottemp=                attfunc("pottemp","potential tempeature (K) at **pref**.")
pottempequisat=         attfunc("pottempequisat","saturation equivalent potential temperature (K) at **pref**.")
pottempequipseudo=      attfunc("pottempequipseudo","pseudo equivalent potential temperature (K) at **pref**.\n\n\
Integrates the pseudo adiabat until the temperature falls\
outside the range of validity of TEOS-10.")
pottempequi=            attfunc("pottempequi","equivalent potential tempeature (K) at **pref**.")

