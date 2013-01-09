========
Tutorial
========
The pyteos_air module is offering a simple interface to the TEOS-10 fortran library.
This tutorial focuses on the ``liq_ice_air`` realm. This realm includes functions
that will work over the whole range of validity of the TEOS-10 wet air component.

The first step is to import the ``liq_ice_air`` realm::

    >>>import pyteos_air.liq_ice_air as liq_ice_air

When this is called the python bindings retrieve the available functions that were interfaced
in the f2py file. These are likely to be expanded in the future. These interfaces carry the i
thermodynamical functions::

    >>>help(liq_ice_air)
   PACKAGE CONTENTS
   g (package)
   g_ref (package)
   h (package)
   il (package)
   sat (package)

For example ``liq_ice_air.g`` carries all the functions in the Gibbs formulation of the wet air::

    >>>help(liq_ice_air.g)
    FUNCTIONS
    cond_entropy(A, T, p)
        :param A: dry air massfraction (kg/kg)
        :type A: np.array.
        :param T: absolute temperature (K)
        :type T: np.array.
        :param p: total pressure (Pa)
        :type p: np.array.
        :returns: wet air entropy with all the moisture in condensed phase (J/K)
        
        This is for wet air with ice and liquid.
    
    entropy(A, T, p)
        :param A: dry air massfraction (kg/kg)
        :type A: np.array.
        :param T: absolute temperature (K)
        :type T: np.array.
        :param p: total pressure (Pa)
        :type p: np.array.
        :returns: wet air entropy (J/K).
        
        At the freezing temperature, this function assumes that half of the condensate is in liquid phase
        and that the other half is in ice phase.
        This is for wet air with ice and liquid.
    
    temperatureequi(A, T, p)
        :param A: dry air massfraction (kg/kg)
        :type A: np.array.
        :param T: absolute temperature (K)
        :type T: np.array.
        :param p: total pressure (Pa)
        :type p: np.array.
        :returns: equivalent temperature (K)

Let's consider the ``liq_ice_air.g.entropy`` function::

    >>>help(liq_ice_air.g.entropy)
    entropy(A, T, p)
    :param A: dry air massfraction (kg/kg)
    :type A: np.array.
    :param T: absolute temperature (K)
    :type T: np.array.
    :param p: total pressure (Pa)
    :type p: np.array.
    :returns: wet air entropy (J/K).
    
    At the freezing temperature, this function assumes that half of the condensate is in liquid phase 
    and that the other half is in ice phase.
    This is for wet air with ice and liquid.

    >>>liq_ice_air.g.entropy(0.9,300,8.5e4)
    array(378.8707010556342)

It returns a numpy array. The pyteos_air functions should handle broadcasting and will accept `np.array`
that have up to three dimensions. These arrays are passed entirely to the Fortran code and there should be
no performance slowdown since the loop is performed in Fortran.

All the other functions can be used similarly.
