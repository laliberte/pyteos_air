from netCDF4 import Dataset
import pyteos_air.liq_ice_air as liq_ice_air
import numpy as np
import pickle
import scipy.interpolate as interp
import multiprocessing as mp

def ccreate_thermo(args):
    #LOAD THE DATA:
    thermo = pickle.load(args.in_thermodynamic_file)
    data = Dataset(args.in_netcdf_file)

    #CREATE THE OUTPUT FILE
    output = Dataset(args.out_netcdf_file,'w',format='NETCDF4',zlib=args.zlib)
    output = replicate_netcdf_file(output,data)

    compute_g(args,thermo,data,output)
    #compute_h(args,thermo,data,output)
    data.close()
    output.close()

def compute_g(args,thermo,data,output):
    #This function computes the functions from interpolated polynomials
    #for input types g, g_ref and sat
    #
    #It is assumed that the input data is using CMIP5 naming convention

    #Pressure and temperature must be in input file:
    var_list=['ta','pa']
    if not set(data.variables.keys()).issuperset(var_list):
        raise IOError('Input file should contain ta and pa')

    #Check that the input coordinates conform:
    shape_list=[]
    for var in var_list:
        shape_list.append(data.variables[var].shape)
    
    if len(set(shape_list))>1:
        raise IOError('Inputs should have the exact same dimensions (this code does not broadcast).')

    #Use temperature as a model variable and create the necessary dimensions for ouput:
    for dims in data.variables['ta'].dimensions:
        if dims not in output.dimension.keys():
            output.createDimension(dims,len(data.dimensions[dims]))
            dim_var = output.createVariable(dims,'d',(dims,))
            dim_var[:] = data.variables[dims][:]
            output = replicate_netcdf_var(output,data,dims)

    #Determine the output variables:
    out_var_list=thermo.keys()
    for var in out_var_list:
        if not thermo[var]._input_type in ['sat','g','g_ref']:
            out_var_list.remove(var)

    if not 'massfraction_air' in out_var_list:
        raise IOError('Variable massfraction_air should be in input thermodynamical file')

    #Create the output variables. Do not output the saturation massfraction
    out_var_list.remove('massfraction_air')
    for var in out_var_list:
        output.createVariable(var,'d',tuple(data.variables['ta'].dimensions))
        output = replicate_netcdf_var_diff(output,data,'ta',var)
        if args.exact>0:
            #If exact was requested (for accuracy testing purposes), create the variables
            output.createVariable(var+'_exact','d',tuple(data.variables['ta'].dimensions))
            output = replicate_netcdf_var_diff(output,data,'ta',var+'_exact')

    #Output the relative humidity if hus is in the inpute file:
    if 'hus' in data.variables.keys():
        output.createVariable('rh_wmo','d',tuple(data.variables['ta'].dimensions))
        output = replicate_netcdf_var_diff(output,data,'ta','rh_wmo')
    
    time_length=len(data.dimensions['time']) 
    for t_id in range(0,time_length):
        T =data.variables['ta'][t_id,...]
        p =data.variables['pa'][t_id,...]

        if 'hus' in data.variables.keys():
            A  =1.0-data.variables['hus'][t_id,...]
            massfraction_air=thermo['massfraction_air'](T,p)
            hur=(1.0 / A - 1.0) / (1.0 / massfraction_air - 1.0)
        elif 'hur' in data.variables.keys():
            hur=data.variables['hur'][t_id,...]
        else:
            raise IOError('Input file should contain hus or hur')

        #First find the relative humidity:
        output.variables['rh_wmo'][t_id,...]=hur
        for var in out_var_list:
            output.variables[var][t_id,...]=thermo[var](hur,T,p)
            if args.exact>0:
                pool=mp.Pool(processes=args.exact)
                output.variables[var+'_exact'][t_id,...]=mp_vec_masked(getattr(getattr(liq_ice_air,thermo[var]._input_type),var),(A,T,p,1e5*np.ones_like(A)),pool=pool)
                pool.close()
        output.sync()
    return

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

def replicate_netcdf_var_diff(output,data,var,var_out):
    for att in data.variables[var].ncattrs():
        if att[0]!='_':
            att_val=getattr(data.variables[var],att)
            if 'encode' in dir(att_val):
                att_val=att_val.encode('ascii','replace')
            setattr(output.variables[var_out],att,att_val)
    return output

def mp_vec_masked(func,args,pool=None):
    fill_value=1e90
    #This function simplifies the use of the multiprocessing toolbox.
    if pool:
        args=np.broadcast_arrays(*args)
        num_procs=len(pool._pool)
        dims_frac_proc=np.ma.array(np.array(args[0].shape)/num_procs)
        dim_index=np.argmin(np.ma.masked_where(dims_frac_proc<1.0,dims_frac_proc))
        in_shape=list(args[0].shape)
        in_shape[dim_index]=1
        iter_list=[[func for x in range(0,args[0].shape[dim_index])]]
        for in_arr in args:
            iter_list.append(np.split(np.ma.filled(in_arr,fill_value),args[0].shape[dim_index],axis=dim_index))
        out_var = np.concatenate(map(lambda x: np.reshape(x,in_shape),
                                    pool.map(tuple_function,zip(*iter_list))),axis=dim_index)
    else:
        out_var = func(*args)
    return np.ma.masked_where(abs(out_var)>=fill_value,out_var)

def tuple_function(args):
    return args[0](*args[1:])

def main():
    import sys
    import argparse 
    import shutil
    import textwrap

    #Option parser
    description=textwrap.dedent('''\
    This script computes the pyteos_air functions stored in a pickle file created by
    create_interpolants.

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

    parser.add_argument('in_thermodynamic_file',type=argparse.FileType('r'),
                         help='Input interpolants file on a rh_wmo, ta, pa grid')
    parser.add_argument('in_netcdf_file',
                         help='Input netCDF file with variables ta, pa, hus or ta, pa, hur. If both sets are available, it uses ta, pa, hus.')
    parser.add_argument('out_netcdf_file',
                         help='Output netCDF file with variables pottemp, pottempequi, rh_wmo')
    parser.add_argument('-z','--zlib',
                         default=False, action='store_true',
                         help='Compress the output using compressed netCDF4.')

    parser.add_argument('-e','--exact',type=int,
                         default=0,
                         help='If set to a value EXACT larger than 0, computes the exact values using liq_ice_air (can be slow) using EAXCT processors.')
    
    args = parser.parse_args()

    create_thermo(args)

if __name__ == "__main__":
    main()
