from netCDF4 import Dataset
import pyteos_air.liq_ice_air as liq_ice_air
import numpy as np
import pickle
import scipy.interpolate as interp
import multiprocessing as mp
import os
import sys
import psutil

#Define the correspondence between pyteos variable names and CMIP5 naming convention:
CMIP5_pyteos_equivalence={'A':'hus',
                          'T':'ta',
                          'p':'pa',
                          'rh_wmo':'hur',
                          'eta':'eta'}
CMIP5_conversions={'A':lambda x: 1.0-x,
                   'T':lambda x: x,
                   'p':lambda x: x,
                   'rh_wmo':lambda x: x,
                   'eta':lambda x: x}

valid_params=[ ('T','p'),
               ('A','T','p'),
               ('rh_wmo','T','p'),
               ('A','eta','p')]

def print_memory_usage(process):
    nr_mbytes = process.get_memory_info()[0] / 1048576.0
    sys.stdout.write("{}\n".format(nr_mbytes))
    sys.stdout.flush()

def create_thermo(args):


    #LOAD THE DATA:
    thermo = pickle.load(args.in_thermodynamic_file)
    data = Dataset(args.in_netcdf_file)

    #CREATE THE OUTPUT FILE
    output = Dataset(args.out_netcdf_file,'w',format='NETCDF4',zlib=args.zlib)
    output = replicate_netcdf_file(output,data)
    fill_value=1e20

    #Determine the output variables:
    out_var_list=thermo.keys()

    #FIRST PASS:
    #Find the available parameters sets:
    available_params=[params for params in valid_params if params_in_data(data,params)]
    #Transfer each of the variables in them to the output file:
    output=transfer_variables(args,data,output,available_params,fill_value)
    data.close()

    variable_list=[]

    while not set(output.variables.keys()).issubset(variable_list):
        variable_list=output.variables.keys()
        available_var_list=[ var for var in out_var_list if function_params(thermo[var])[:3] in available_params]
        for var in available_var_list:
            output = create_output(args,output,thermo[var],fill_value)

        #CREATE rh_wmo IF massfraction_air IS AVAILABLE:
        if params_in_output(output,('A','massfraction_air')) and 'rh_wmo' not in output.variables.keys():
            #rh_wmo_function=np.vectorize(lambda A, massfraction_air: np.where(massfraction_air>0.0,(1.0 / A - 1.0) / (1.0 / massfraction_air - 1.0),0.0))
            rh_wmo_function=(lambda A, massfraction_air: (1.0 - A) / (1.0 - massfraction_air) * (massfraction_air / A) )
            rh_wmo_function.__name__='rh_wmo'
            rh_wmo_function.__doc__='rh_wmo(A,massfraction_air)'
            output = create_output(args,output,rh_wmo_function,fill_value)
            del rh_wmo_function

        if params_in_output(output,('rh_wmo','massfraction_air')) and 'A' not in output.variables.keys():
            #rh_wmo_function=np.vectorize(lambda rh_wmo, massfraction_air: np.where(massfraction_air>0.0,1.0 / (1.0 + rh_wmo * (1.0 / massfraction_air - 1.0)),1.0))
            rh_wmo_function=(lambda rh_wmo, massfraction_air: massfraction_air / (massfraction_air + rh_wmo * (1.0 - massfraction_air)))
            rh_wmo_function.__name__='A'
            rh_wmo_function.__doc__='A(rh_wmo,massfraction_air)'
            output = create_output(args,output,rh_wmo_function,fill_value)
            del rh_wmo_function

        #SECOND PASS:
        #Find the available parameters sets:
        available_params=[params for params in valid_params if params_in_output(output,params)]

        available_var_list=[ var for var in out_var_list if function_params(thermo[var])[:3] in available_params]
        for var in available_var_list:
            output = create_output(args,output,thermo[var],fill_value)

    output.close()

def params_in_data(data,params):
    return set([CMIP5_pyteos_equivalence[var] for var in params]).issubset(data.variables.keys())

def params_in_output(data,params):
    return set(params).issubset(data.variables.keys())

def transfer_variables(args,data,output,available_params,fill_value):
    for params in available_params:
        for var in params:
            for dims in data.variables[CMIP5_pyteos_equivalence[var]].dimensions:
                if dims not in output.dimensions.keys():
                    output.createDimension(dims,len(data.dimensions[dims]))
                    dim_var = output.createVariable(dims,'d',(dims,))
                    dim_var[:] = data.variables[dims][:]
                    output = replicate_netcdf_var(output,data,dims)
                    output.sync()
            if var not in output.variables.keys():
                coord_var=output.createVariable(var,'f',tuple(data.variables[CMIP5_pyteos_equivalence[var]].dimensions),fill_value=fill_value,zlib=args.zlib)
                output = replicate_netcdf_var_diff(output,data,CMIP5_pyteos_equivalence[var],var)
                coord_var[:] = CMIP5_conversions[var](data.variables[CMIP5_pyteos_equivalence[var]][:])
                output.sync()
    return output

def create_output(args,output,func,fill_value):
    if not func.__name__ in output.variables.keys():
        params_list=function_params(func)[:3]
        output.createVariable(func.__name__,'f',tuple(output.variables[params_list[0]].dimensions),fill_value=fill_value,zlib=args.zlib)
        
        time_length=len(output.dimensions['time']) 
        for t_id in range(0,time_length):
            coordinates=[]
            for var in params_list:
                coordinates.append(output.variables[var][t_id,...])

            output.variables[func.__name__][t_id,...]=np.ma.filled(mp_vec_masked(func,coordinates),fill_value=fill_value)
            del coordinates
            output.sync()
    #print_memory_usage(args.process)
    return output

#    if args.exact>0:
#        output.createVariable(func.__name__+'_exact','d',tuple(output.variables[params_list[0]].dimensions),fill_value=fill_value)
#        if args.exact>0:
#            pool=mp.Pool(processes=args.exact)
#            if func._input_type=='g_ref':
#                coordinates.append(1e5*np.ones_like(coordinates[0]))
#            output.variables[func.__name__+'_exact'][t_id,...]=np.ma.filled(
#                                                        mp_vec_masked(getattr(
#                                                                         getattr(
#                                                                            liq_ice_air,func._input_type),func.__name__
#                                                                                ),
#                                                                         tuple(coordinates),
#                                                                         pool=pool
#                                                                     ),fill_value=fill_value
#                                                                  )
#            pool.close()

def function_params(interp_func):
    return tuple(interp_func.__doc__.splitlines()[0].replace(interp_func.__name__,'').replace('(','').replace(')','').split(','))

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
    fill_value=1e20
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
    #out_var = np.ma.masked_where(np.isnan(out_var),out_var)
    #return np.ma.masked_where(np.abs(out_var)>=fill_value,out_var)
    return np.ma.fix_invalid(out_var,fill_value=fill_value)

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

    from the air temperature (ta), air pressure (pa), and the Specific 
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
                         help='If set to a value EXACT larger than 0, computes the exact values using liq_ice_air (can be extremely slow) using EXACT processors.')
    
    args = parser.parse_args()

    args.process = psutil.Process(os.getpid())
    create_thermo(args)

if __name__ == "__main__":
    main()
