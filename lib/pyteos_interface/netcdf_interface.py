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

    data = Dataset(args.in_netcdf_file)

    #CREATE THE OUTPUT FILE
    output = Dataset(args.out_netcdf_file,'w',format='NETCDF4',zlib=args.zlib)
    output = replicate_netcdf_file(output,data)

    var_list=['ta','pa','hus']
    if not set(data.variables.keys()).issuperset(var_list):
        raise IOError('Input file should contain ta, pa, hus')

    shape_list=[]
    for var in var_list:
        shape_list.append(data.variables[var].shape)
    
    if len(set(shape_list))>1:
        raise IOError('Variables ta, pa and hus should have the exact same dimensions (this code does not broadcast).')

    for dims in data.variables['ta'].dimensions:
        output.createDimension(dims,len(data.dimensions[dims]))
        dim_var = output.createVariable(dims,'d',(dims,))
        dim_var[:] = data.variables[dims][:]
        output = replicate_netcdf_var(output,data,dims)

    out_var_list=thermo.variables.keys()
    if not set(['rh_wmo','a_sat']).issubset(out_var_list):
        raise IOError('Variable rh_wmo and a_sat should be in input thermodynamical file')
    for var in out_var_list:
        thermo_data[var]=thermo.variables[var][:]
        if var in dim_list.keys():
            thermo_data[var]=np.reshape(thermo_data[var],dim_list[var])
        else:
            output.createVariable(var,'f',data.variables['ta'].dimensions)
    #Deal with rh_wmo explicitly:
    out_var_list.remove('rh_wmo')
    out_var_list.remove('a_sat')
    #Remove dimensions:
    out_var_list=[ x for x in out_var_list if (x not in dim_list)]
        
    interp={}
    interp['a_sat']=interp_linear(thermo_data['a_sat'],thermo_data['ta'],thermo_data['pa'])
    for var in out_var_list:
        print var
        interp['var']=interp_linear(thermo_data[var],
             thermo_data['hur'],
             thermo_data['ta'],
             thermo_data['pa']
             )

    time_length=len(data.dimensions['time']) 
    for t_id in range(0,time_length):
        print(t_id)
        a  =1.0-data.variables['hus'][t_id,...]
        ta =data.variables['ta'][t_id,...]
        pa =data.variables['pa'][t_id,...]

        ax_ind=np.argmax(ta.shape)
        #First find the relative humidity:
        a_sat=np.reshape(a_sat['a_sat'](coord_trans(ta,pa)),ta.shape)
        hur=(1.0 / a - 1.0) / (1.0 / a_sat - 1.0)
        output.variables['rh_wmo'][t_id,...]=hur
        for var in out_var_list:
            data.variables[var][t_id,...]=np.reshape(
                                                     interp[var](coors_trans(hur,ta,pa)),
                                                     ta.shape
                                                     )
        output.sync()

    output.close()

def interp_linear(nd_array,*args):
    nd_interp=interp.LinearNDInterpolator(
                            coord_trans(*args),
                            np.ravel(nd_array))
    return nd_interp

def coord_trans(*arr):
    return np.vstack(map(lambda x: np.ravel(x),np.broadcast_arrays(*arr))).T

def replicate_netcdf_file(output,data):
	for att in data.ncattrs():
            att_val=getattr(data,att)
            if 'encode' in dir(att_val):
                att_val=att_val.encode('ascii','replace')
	    setattr(output,att,att_val)
        #output.history+='\n' 
        #output.history+=dt.datetime.now().strftime('%Y-%m-%d %H:%M') #Add time
        #output.history+=' joint_distribution.py'
	return output

def replicate_netcdf_var(output,data,var):
	for att in data.variables[var].ncattrs():
	    if att[0]!='_':
                att_val=getattr(data.variables[var],att)
                if 'encode' in dir(att_val):
                    att_val=att_val.encode('ascii','replace')
	        setattr(output.variables[var],att,att_val)
	return output

def replicate_netcdf_var2(output,data,var,var_out):
	for att in data.variables[var].ncattrs():
	    if att[0]!='_':
                att_val=getattr(data.variables[var],att)
                if 'encode' in dir(att_val):
                    att_val=att_val.encode('ascii','replace')
	        setattr(output.variables[var_out],att,att_val)
	return output

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
    parser.add_argument('in_netcdf_file',
                         help='Input netCDF file with variables ta, pa, hus')
    parser.add_argument('out_netcdf_file',
                         help='Output netCDF file with variables pottemp, pottempequi, rh_wmo')
    parser.add_argument('-z','--zlib',
                         default=False, action='store_true',
                         help='Compress the output using compressed netCDF4.')
    
    args = parser.parse_args()

    compute_both_theta(args)

