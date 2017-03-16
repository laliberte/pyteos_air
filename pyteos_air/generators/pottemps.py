#Standard packages:
import pickle
import numpy as np
import click

#Specific packages:
from .. import liq_ice_air
from ..pyteos_interface import create_gridded_data
from ..create_interpolants import Interpolated_data

@click.command()
@click.argument('dest_dir')
@click.option('--num_procs', type=int,
              default=1,
              help="Number of processes to use. Default: 1")
@click.option('--test', is_flag=True,
              help="Simple testing. Create a much less precise table")
def compute_pottemps(dest_dir, num_procs, test):
    """
    \b
    This script computes two data files:
        DEST_DIR/massfraction_air_pottempequi_g_common_pressures.dat  
        DEST_DIR/massfraction_air_pottemp_g_ref_common_pressures.dat

    that can be used with pyteos_nc to compute moist air entropy and the dry (virtual) potential temperature
    from CF-compliants netcdf files.
    The netcdf files should contain variables ta (air temperature), hus (specific humidity) and air pressure (pa).
    """
    name = 'common_pressures'
    if test:
        thermo_axes = {'rh_wmo': np.linspace(0, 2.0, 11),
                       'T': np.linspace(193, 330, 23),
                       'p': np.array([250.0,500.0,850.0,1000.0])*1e2}
    else:
        thermo_axes = {'rh_wmo': np.linspace(0, 2.0, 101),
                       'T': np.linspace(193, 330, 138),
                       'p': np.array([250.0,300.0,400.0,500.0,600.0,700.0,850.0,925.0,1000.0])*1e2}
    compute_pottempequi(num_procs, thermo_axes, name, dest_dir)
    thermo_axes['pref'] = 1e5
    compute_pottempliquid(num_procs, thermo_axes, name, dest_dir)


def compute_pottempequi(num_procs, thermo_axes, name, dest_dir):
    func_list = {}
    functions = ['pottempequi']
    realm = 'g'
    for function in functions:
        if function == 'pottempequi':
            thermo_data = create_gridded_data(liq_ice_air, realm, 'entropy', thermo_axes,
                                              num_procs=num_procs)
            thermo_data = 273.15 * np.exp(thermo_data/1003.0)
        else:
            thermo_data = create_gridded_data(liq_ice_air, realm, function, thermo_axes,
                                              num_procs=num_procs)
        func_list[function] = Interpolated_data(liq_ice_air, realm, function, thermo_axes, thermo_data)

    func_list['massfraction_air'] = compute_massfraction_air(num_procs, thermo_axes)

    with open(dest_dir+'/'+'_'.join(func_list.keys()) + '_' + realm + '_' + name + '.dat', 'w') as file:
        pickle.dump(func_list,file)


def compute_massfraction_air(num_procs, thermo_axes):
    thermo_axes_massfrac = {'T': thermo_axes['T'],
                            'p': thermo_axes['p']}
    thermo_data = create_gridded_data(liq_ice_air, 'sat', 'massfraction_air',
                                      thermo_axes_massfrac, num_procs=num_procs)
    return Interpolated_data(liq_ice_air, 'sat', 'massfraction_air', thermo_axes_massfrac,
                             thermo_data)


def compute_pottempliquid(num_procs, thermo_axes, name, dest_dir):
    func_list = {}
    functions=['pottemp']
    realm='g_ref'

    for function in functions:
        thermo_data = create_gridded_data(liq_ice_air,realm,function,thermo_axes,num_procs=num_procs)
        func_list[function] = Interpolated_data(liq_ice_air,realm,function,thermo_axes,thermo_data)

    func_list['massfraction_air'] = compute_massfraction_air(num_procs, thermo_axes)

    with open(dest_dir+'/'+'_'.join(func_list.keys()) + '_' + realm + '_' + name + '.dat', 'w') as file:
        pickle.dump(func_list,file)


if __name__=='__main__':
    compute_thermo()
