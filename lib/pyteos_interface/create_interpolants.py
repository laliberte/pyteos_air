import numpy as np
import scipy.interpolate as interp

from create_gridded_data import create_gridded_data

class Interpolated_data:
    def __init__(self,realm,input_type,func_name,thermo_axes,num_procs=1):
        #self._realm=realm
        self._input_type=input_type
        self.__name__=getattr(getattr(getattr(realm,input_type),func_name),'__name__')

        params_list=function_parameters(getattr(getattr(realm,input_type),func_name))
        if input_type in ['g','g_ref'] and 'rh_wmo' in thermo_axes.keys():
            params_list[0]='rh_wmo'
        docstring=self.__name__+'('+','.join(params_list)+')'+getattr(getattr(getattr(realm,input_type),func_name),'__doc__')
        docstring+='\nSpline interpolation from gridded data with grid:\n'

        thermo_data=create_gridded_data(realm,input_type,func_name,thermo_axes,num_procs=num_procs)

        if input_type=='sat':
            docstring+='T: '+str(thermo_axes['T'])+'\n'
            docstring+='p: '+str(thermo_axes['p'])+'\n'

            self._interpolants=interp_linear((thermo_axes['T'],thermo_axes['p'],np.squeeze(thermo_data)))
        elif input_type in ['h']:
            self._p = thermo_axes['p']
            num_p=len(self._p)

            docstring+='A: '+str(thermo_axes['A'])+'\n'
            docstring+='eta: '+str(thermo_axes['eta'])+'\n'
            docstring+='p: '+str(thermo_axes['p'])+'\n'
            docstring+='The interpolation is a spline on constant p and linear along p.'

            self._interpolants=interp_linear_3D((thermo_axes['A'],thermo_axes['eta'],thermo_axes['p'],thermo_data))
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


            self._interpolants=interp_linear_3D((thermo_axes[moisture_var],thermo_axes['T'],thermo_axes['p'],thermo_data))
        self.__doc__=docstring

    def __call__(self,*args):
        if self._input_type=='sat':
            T, p = np.broadcast_arrays(*np.atleast_3d(*args))
            return np.reshape(self._interpolants.ev(np.ravel(T),np.ravel(p)),T.shape)
        elif self._input_type in ['h','g','g_ref']:
            return self._interpolants(np.concatenate([a[...,np.newaxis] for a in np.broadcast_arrays(*np.atleast_3d(*args))],axis=-1))
        else:
            return
        
def interp_linear(args):
    nd_interp=interp.RectBivariateSpline(*args,kx=1,ky=1)
    return nd_interp

def interp_linear_3D(args):
    nd_interp=interp.RegularGridInterpolator(args[:-1],args[-1],bounds_error=False)
    return nd_interp

def function_parameters(func):
    func_desc=getattr(func,'func_doc')
    return [ x.split(':')[1][6:] for x in func_desc.split('\n') if (len(x)>=6 and x[:6]==':param')]

def pickle_interpolants(realm,input_type,functions_list,thermo_axes,file_name,num_procs=1):
    """
    :param realm: Realm (e.g. 'lid_ice_air')
    :type realm: str.
    :param input_type: Type of arguments (e.g. 'g')
    :type T: str.
    :param functions_list: List of function (e.g. ['density','enthalpy'])
    :type functions_list: list.
    :param thermo_axes: dictionary of axes to interpolate to (e.g. {'A':[0.5,1.0],'T':[250,300],'p':[1e4,1e5]})
    :type thermo_axes: dict.
    :param file_name: Name of file to pickle to
    :type file_name: str.

    This function creates a pickle file that can be used with pyteos_nc. This greatly speeds up computation.
    """

    for function in functions:
        func_list[function]=Interpolated_data(liq_ice_air,realm,function,thermo_axes,num_procs=num_procs)
    file=open(file_name,'w')
    pickle.dump(func_list,file)

