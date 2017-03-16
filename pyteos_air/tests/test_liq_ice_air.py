"""
    Test the table generator.
"""

import numpy as np
from contextlib import closing

from pyteos_air import liq_ice_air


def test_g_entropy():
    # Below freezing:
    assert liq_ice_air.g.entropy(0.98, 270.0, 1e5) == 167.3697059091392
    # Above freezing:
    assert liq_ice_air.g.entropy(0.98, 300.0, 1e5) == 276.53730016658955
    # Low pressure, must be dry:
    assert liq_ice_air.g.entropy(1.0, 300.0, 5e3) == 958.6649347968464


def test_g_ref_pottemp():
    # Below freezing:
    assert liq_ice_air.g_ref.pottemp(0.98, 270.0, 5e4, 1e5) == 328.77732775760836
    # Above freezing:
    assert liq_ice_air.g_ref.pottemp(0.98, 300.0, 5e4, 1e5) == 365.2154567609769
    # Low pressure, must be dry:
    assert liq_ice_air.g_ref.pottemp(1.0, 300.0, 5e3, 1e5) == 692.815889987961


def test_g_ref_pottempequiapprox():
    # Below freezing:
    assert liq_ice_air.g_ref.pottempequiapprox(0.98, 270.0, 5e4, 1e5) == 395.17024568743096
    # Above freezing:
    assert liq_ice_air.g_ref.pottempequiapprox(0.98, 300.0, 5e4, 1e5) == 440.05599314639227
    # Low pressure, must be dry:
    assert liq_ice_air.g_ref.pottempequiapprox(1.0, 300.0, 5e3, 1e5) == 710.3932123977726
    # Above freezing at surface:
    assert liq_ice_air.g_ref.pottempequiapprox(0.98, 300.0, 1e5, 1e5) == 359.86582890653585
