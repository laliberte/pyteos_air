import numpy as np

def function_parameters(func):
    func_desc=getattr(func,'func_doc')
    return [ x.split(':')[1][6:] for x in func_desc.split('\n') if (len(x)>=6 and x[:6]==':param')]

def create_gridded_data(realm,input_type,func,thermo_axes,num_procs=1):
    params_list=function_parameters(getattr(getattr(realm,input_type),func))
    if input_type not in ['g','g_ref','h','sat']:
        raise NotImplementedError('Input type '+input_type+' is not implemented.')

    #Check that the number of function parameters equals the number of thermo_axes inputs:
    if len(params_list) != len(thermo_axes.keys()):
        raise NameError('class Interpolated data: the input thermodynamics axes must have the same number of axes as the functions inputs.')

    #Check that all the functions parameters are given:
    if 'A' in params_list:
            fixed_params_list=[ x if x is not 'A' else 'rh_wmo' for x in params_list]
    if not set(params_list).issubset(thermo_axes.keys()) and not set(fixed_params_list).issubset(thermo_axes.keys()):
        raise NameError('class Interpolated data: the input thermodynamics axes must contain all of the functions inputs.')

    if input_type=='g_ref' and not np.isscalar(thermo_axes['pref']):
        raise IndexError('class Interpolated data: if the function has more than three dimensions pref must be a scalar')

    #Start the multiprocessing pool:
    if num_procs>1:
        import multiprocessing as mproc
        pool=mproc.Pool(processes=int(num_procs))
    else:
        pool=None

    if input_type in ['g','g_ref','sat']:
        T =np.reshape(thermo_axes['T'],[1,len(thermo_axes['T']),1])
        p =np.reshape(thermo_axes['p'],[1,1,len(thermo_axes['p'])])

        if input_type in ['g','g_ref']:
            if 'rh_wmo' in thermo_axes.keys():
                rh_wmo = np.reshape(thermo_axes['rh_wmo'],[len(thermo_axes['rh_wmo']),1,1])
                a_sat=np.reshape(mp_vec_masked(getattr(getattr(realm,'sat'),'massfraction_air'),(T,p)),(1,len(thermo_axes['T']),len(thermo_axes['p'])))
                A = 1.0 / (1.0 + rh_wmo * (1.0 / a_sat - 1.0))
                #Fix problems close to saturation:
                A = np.where(rh_wmo>=1.0,A+1e-10,A)
            else:
                A = np.reshape(thermo_axes['A'],[len(thermo_axes['A']),1,1])
            #If temperature is less than 193.0, make it dry:
            A=np.where(T<193.0,1.0,A)

    elif input_type in ['h']:
        eta =np.reshape(thermo_axes['eta'],[1,len(thermo_axes['eta']),1])
        p =np.reshape(thermo_axes['p'],[1,1,len(thermo_axes['p'])])
        A = np.reshape(thermo_axes['A'],[len(thermo_axes['A']),1,1])
        
    if input_type=='g_ref':
        result=mp_vec_masked(getattr(getattr(realm,input_type),func),(A,T,p,thermo_axes['pref']),pool=pool)
    elif input_type=='g':
        result=mp_vec_masked(getattr(getattr(realm,input_type),func),(A,T,p),pool=pool)
    elif input_type=='sat':
        result=mp_vec_masked(getattr(getattr(realm,input_type),func),(T,p),pool=pool)
    elif input_type=='h':
        result=mp_vec_masked(getattr(getattr(realm,input_type),func),(A,eta,p),pool=pool)

    if num_procs>1:
        pool.close()

    return result

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

def pickle_function(realm,input_type,func,thermo_axes,file_to_pickle,num_procs=1):
    import pickle
    import pyteos_air.liq_ice_air as liq_ice_air
    func_instance=Interpolated_data(realm,input_type,func,thermo_axes,num_procs=num_procs)

    file_id=open(file_to_pickle,'w')
    pickle.dump(func_instance,file_id)
    file_id.close
