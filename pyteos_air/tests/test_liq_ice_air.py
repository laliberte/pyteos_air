"""
    Test the table generator.
"""

import numpy as np
from contextlib import closing

from pyteos_air import liq_ice_air


def test_g_entropy():
    # From src/SIA_library/Air_3b.F90
    # assert liq_ice_air.g.entropy(0.9, 300.0, 1e5) == 911.170080461
    # From src/SIA_library/Ice_Air_4b.F90
    assert liq_ice_air.g.entropy(0.5, 270.0, 1e5) == -610.264515318
    # From src/SIA_library/Liq_Air_4b.F90
    assert liq_ice_air.g.entropy(0.5, 300.0, 1e5) == 343.783393872
