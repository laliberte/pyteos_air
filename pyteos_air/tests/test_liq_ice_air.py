"""
    Test the table generator.
"""

import numpy as np
from contextlib import closing

from pyteos_air import liq_ice_air


def test_g_entropy():
    # From src/SIA_library/Air_3b.F90
    assert liq_ice_air.g.entropy(0.9, 300.0, 1e5) == 911.170080461
