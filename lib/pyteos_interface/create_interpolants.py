import numpy as np
import scipy.interpolate as interp
import matplotlib.pyplot as plt

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

            self._interpolants=map(
                            interp_linear,
                            zip(
                                [thermo_axes['A'] for x in range(0,num_p)],
                                [thermo_axes['eta'] for x in range(0,num_p)],
                                np.dsplit(thermo_data,num_p)
                            )
                            )
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
            T, p = np.broadcast_arrays(*np.atleast_3d(*args))
            return np.reshape(self._interpolants.ev(np.ravel(T),np.ravel(p)),T.shape)
        elif self._input_type in ['h','g','g_ref']:
            A, T, p =np.broadcast_arrays(*np.atleast_3d(*args))
            #Find the pressure bin:
            p_ind_sup=np.apply_along_axis(np.digitize,0,p,self._p)
            p_ind_low=p_ind_sup-1
            p_ind_sup=np.minimum(p_ind_sup,len(self._p)-1)
            p_ind_low=np.maximum(p_ind_low,0)
        
            pval=np.vectorize(lambda x: self._p[x])
            p_sup=pval(p_ind_sup)
            p_low=pval(p_ind_low)
            dp=p_sup-p_low
        
            interpolants=np.vectorize(lambda x,y,z: self._interpolants[z].ev(x,y))
            return interpolants(A,T,p_ind_sup)*(p-p_low)/dp + interpolants(A,T,p_ind_low)*(p_sup-p)/dp
        else:
            return
        
def interp_linear(args):
    nd_interp=interp.RectBivariateSpline(*args,kx=1,ky=1)
    return nd_interp

def function_parameters(func):
    func_desc=getattr(func,'func_doc')
    return [ x.split(':')[1][6:] for x in func_desc.split('\n') if (len(x)>=6 and x[:6]==':param')]

