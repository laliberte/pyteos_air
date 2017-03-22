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
    assert liq_ice_air.g_ref.pottempequiapprox(0.98, 270.0, 5e4, 1e5) == 393.6800617325973
    # Above freezing:
    assert liq_ice_air.g_ref.pottempequiapprox(0.98, 300.0, 5e4, 1e5) == 438.39654538337976
    # Low pressure, must be dry:
    assert liq_ice_air.g_ref.pottempequiapprox(1.0, 300.0, 5e3, 1e5) == 707.7143250617679
    # Above freezing at surface:
    assert liq_ice_air.g_ref.pottempequiapprox(0.98, 300.0, 1e5, 1e5) == 358.5087776356424
    # Above freezing at surface mostly dry:
    assert liq_ice_air.g_ref.pottempequiapprox(0.999, 300.0, 1e5, 1e5) == 303.1831764842107


def test_g_ref_pottempequi_exact():
    # Below freezing:
    assert liq_ice_air.g_ref.pottempequi(0.98, 270.0, 5e4, 1e5) == 383.8781880195045
    # Above freezing:
    assert liq_ice_air.g_ref.pottempequi(0.98, 300.0, 5e4, 1e5) == 424.2285658682927
    # Low pressure, must be dry:
    assert liq_ice_air.g_ref.pottempequi(1.0, 300.0, 5e3, 1e5) == 692.815889987961
    # Above freezing at surface:
    assert liq_ice_air.g_ref.pottempequi(0.98, 300.0, 1e5, 1e5) == 351.77280280400373
    # Above freezing at surface mostly dry:
    assert liq_ice_air.g_ref.pottempequi(0.999, 300.0, 1e5, 1e5) == 302.37955703334933
