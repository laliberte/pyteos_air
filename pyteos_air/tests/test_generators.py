"""
    Test the table generator.
"""

import numpy as np
from contextlib import closing

from pyteos_air.generators.pottemps import compute_pottemps


#def test_pottemps(tmpdir):
#    dest_dir = tmpdir.mkdir('data')
#    compute_pottemps(str(dest_dir), 1, True) 
#    assert 
