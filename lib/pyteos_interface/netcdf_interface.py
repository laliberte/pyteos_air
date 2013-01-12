from netCDF4 import Dataset
import pyteos_air.liq_ice_air as liq_ice_air
import numpy as np
import pickle
import scipy.interpolate as interp

def compute_both_theta(args):
    #LOAD THE DATA:
    thermo = Dataset(args.in_thermodynamic_file)
    thermo_data=dict()
    dim_list={'hur':[len(thermo.dimensions['hur']),1,1],
                'ta':[1,len(thermo.dimensions['ta']),1],
                'pa':[1,1,len(thermo.dimensions['pa'])]}

    var_list=['ta','pa','hus']

    out_var_list=thermo.variables.keys()
    if not set(['rh_wmo','a_sat']).issubset(out_var_list):
        raise IOError('Variable rh_wmo and a_sat should be in input thermodynamical file')
    for var in out_var_list:
        thermo_data[var]=thermo.variables[var][:]
        if var in dim_list.keys():
            thermo_data[var]=np.reshape(thermo_data[var],dim_list[var])
    #Deal with rh_wmo explicitly:
    out_var_list.remove('rh_wmo')
    out_var_list.remove('a_sat')
    #Remove dimensions:
    out_var_list=[ x for x in out_var_list if (x not in dim_list)]
        
    interp={}
    interp['a_sat']=interp_linear(thermo_data['a_sat'],thermo_data['ta'],thermo_data['pa'])
    pickle.dump(interp,args.out_pickle_file)
    #for var in out_var_list:
    #    print var
    #    interp['var']=interp_linear(thermo_data[var],
    #         thermo_data['hur'],
    #         thermo_data['ta'],
    #         thermo_data['pa']
    #         )
    #pickle.dump(interp,args.out_pickle_file)

def interp_linear(nd_array,*args):
    nd_interp=interp.interpnd(
                            coord_trans(*args),
                            np.ravel(nd_array))
    return nd_interp

def coord_trans(*arr):
    return np.vstack(map(lambda x: np.ravel(x),np.broadcast_arrays(*arr))).T

if __name__ == "__main__":
    import sys
    import argparse 
    import shutil
    import textwrap

    #Option parser
    description=textwrap.dedent('''\
    This script computes:
    1. Liquid water (Dry) Potential Temperature, output as pottemp
    2. Equivalent (Moist) Potential Temperature, output as pottempequi
    3. Relative Humidity using the WMO definition, output as rh_wmo

    from the air temperature (a), air pressure (pa), and the Specific 
    Total Water Content (hus) which can be taken equal to the
    Specific Humidity if the climate model does not allow condensated
    water.
    ''')
    epilog='Frederic Laliberte, Paul Kushner 11/2012'
    epilog=textwrap.dedent(epilog+'\n\nThis script uses the pyteos_air library, based on TEOS-10.')
    version_num='0.1'
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                            description=description,
                            version='%(prog)s '+version_num,
                            epilog=epilog)

    parser.add_argument('in_thermodynamic_file',
                         help='Input netCDF file with variables pottemp,potempequi, rh_wmo on a rh_wmo, ta, pa grid')
    parser.add_argument('out_pickle_file',type=argparse.FileType('w'),
                         help='Output netCDF file with variables pottemp, pottempequi, rh_wmo')
    
    args = parser.parse_args()

    compute_both_theta(args)

