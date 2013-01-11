import numpy as np
import multiprocessing as mproc
import pickle
import scipy.interpolate as interp

fill_value=1e90

class Interpolated_data:
    def __init__(self,realm,input_type,func,thermo_axes,num_procs=1):
        self._realm=realm
        self._input_type=input_type
        self.__name__=getattr(getattr(getattr(realm,input_type),func),'__name__')
        params_list=function_parameters(getattr(getattr(realm,input_type),func))

        docstring=self.__name__+'('+','.join(params_list)+')'+getattr(getattr(getattr(realm,input_type),func),'__doc__')
        docstring+='\nSpline interpolation from gridded data with grid:\n'

        if input_type not in ['g','g_ref','sat']:
            raise NotImplementedError('Input type '+input_type+'is not implemented.')

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

        thermo_data = create_gridded_data(realm,input_type,func,thermo_axes,num_procs)

        if input_type=='sat':
            docstring+='T: '+str(thermo_axes['T'])+'\n'
            docstring+='p: '+str(thermo_axes['p'])+'\n'

            self._interpolants=interp_linear((thermo_axes['T'],thermo_axes['p'],thermo_data))
        elif input_type in ['g','g_ref']:
            self._p = thermo_axes['p']
            num_p=len(self._p)

            if 'rh_wmo' in thermo_axes.keys():
                moisture_var='rh_wmo'
                docstring=docstring.replace(':param A: dry air massfraction (kg/kg)',
                                            ':param rh_wmo: relative humidity, WMO definition, 0 < rh < 1 for unsaturated air')  
                docstring=docstring.replace(':type A:',
                                            ':type rh_wmo')
                docstring+='rh_wmo: '+str(thermo_axes['rh_wmo'])+'\n'
            else:
                moisture_var='A'
                docstring+='A: '+str(thermo_axes['A'])+'\n'

            docstring+='T: '+str(thermo_axes['T'])+'\n'
            docstring+='p: '+str(thermo_axes['p'])+'\n'
            docstring+='The interpolation is a spline on constant p and linear along p.'

            self._interpolants=map(
                            interp_linear,
                            zip(
                                [thermo_axes[moisture_var] for x in range(0,num_p)],
                                [thermo_axes['T'] for x in range(0,num_p)],
                                np.dsplit(thermo_data,num_p)
                            )
                            )
        self.__doc__=docstring

    def __call__(self,*args):
        
        if self._input_type=='sat':
            T, p = np.atleast_3d(*args)
            return self._interpolants(T,p)
        elif self._input_type in ['g','g_ref']:
            rh_wmo, T, p =np.atleast_3d(*args)
            #Find the pressure bin:
            p_ind_sup=np.apply_along_axis(np.digitize,0,p,self._p)
            p_ind_low=p_ind_sup-1
            p_ind_sup=np.minimum(p_ind_sup,len(self._p)-1)
            p_ind_low=np.maximum(p_ind_low,0)

            p_sup=self._p[p_ind_sup]
            p_low=self._p[p_ind_low]
            dp=p_sup-p_low

            interpolants=np.vectorize(lambda x,y,z: self._interpolants[z](x,y))
            return interpolants(rh_wmo,T,p_ind_sup)*(p-p_low)/dp + interpolants(rh_wmo,T,p_ind_low)*(p_sup-p)/dp
        
def interp_linear(args):
    nd_interp=interp.RectBivariateSpline(*args)
    return nd_interp

def function_parameters(func):
    func_desc=getattr(func,'func_doc')
    return [ x.split(':')[1][6:] for x in func_desc.split('\n') if (len(x)>=6 and x[:6]==':param')]

def create_gridded_data(realm,input_type,func,thermo_axes,num_procs=1):
    params_list=function_parameters(getattr(getattr(realm,input_type),func))

    T =np.reshape(thermo_axes['T'],[1,len(thermo_axes['T']),1])
    p =np.reshape(thermo_axes['p'],[1,1,len(thermo_axes['p'])])

    if input_type in ['g','g_ref']:
        if 'rh_wmo' in thermo_axes.keys():
            rh_wmo = np.reshape(thermo_axes['rh_wmo'],[len(thermo_axes['rh_wmo']),1,1])
            a_sat=np.reshape(getattr(getattr(realm,'sat'),'massfraction_air')(T,p),(1,len(thermo_axes['T']),len(thermo_axes['p'])))
            A = 1.0 / (1.0 + rh_wmo * (1.0 / a_sat - 1.0))
            #Fix problems close to saturation:
            A = np.where(rh_wmo>=1.0,A+1e-10,A)
        else:
            A = np.reshape(thermo_axes['A'],[len(thermo_axes['A']),1,1])

    if num_procs>1:
        pool=mproc.Pool(processes=int(num_procs))
    else:
        pool=None

    if input_type=='g_ref':
        result=mp_vec_masked(getattr(getattr(realm,input_type),func),(A,T,p,thermo_axes['pref']),pool=pool)
    elif input_type=='g':
        result=mp_vec_masked(getattr(getattr(realm,input_type),func),(A,T,p),pool=pool)
    elif input_type=='sat':
        result=mp_vec_masked(getattr(getattr(realm,input_type),func),(A,T),pool=pool)

    if num_procs>1:
        pool.close()

    return result

def mp_vec_masked(func,args,pool=None):
    #This function simplifies the use of the multiprocessing toolbox.
    if pool:
        num_procs=len(pool._pool)
        dims_frac_proc=np.ma.array(np.array(args[0].shape)/num_procs)
        dim_index=np.argmin(np.ma.masked_where(dims_frac_proc<1.0,dims_frac_proc))
        iter_list=[[func for x in range(0,args[0].shape[dim_index])]]
        for in_arr in args:
            iter_list.append([np.ma.filled(in_arr,fill_value) for x in range(0,in_arr.shape[dim_index])])
        out_var = np.concatenate(pool.map(tuple_function,zip(*iter_list)),axis=dim_index)
    else:
        out_var = func(*args)
    return np.ma.masked_where(abs(out_var)>=fill_value,out_var)

def tuple_function(args):
    return args[0](*args[1:])
