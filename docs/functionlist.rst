=============
Function List
=============
Realm ``liq_ice_air``
---------------------
This realm takes into account phase changes. For the Gibbs formulations, it assumes that at
the freezing point half of the condensate is ice and the other half is liquid water. This
is an arbitray constraint that is not necessary in the enthalpy formulation. The enthalpy formulation
is thus fully determined by physical constraints but not the Gibbs formulation.

The Gibbs formulation ``g``
^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. automodule:: pyteos_air.liq_ice_air.g
    :members:

The Gibbs formulation with reference values ``g_ref``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. automodule:: pyteos_air.liq_ice_air.g_ref
    :members:

The enthalpy formulation ``h``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. automodule:: pyteos_air.liq_ice_air.h
    :members:

The saturation values ``sat``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. automodule:: pyteos_air.liq_ice_air.sat
    :members:

Isentropic levels ``il``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. automodule:: pyteos_air.liq_ice_air.il
    :members:

