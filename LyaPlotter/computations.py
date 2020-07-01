# %% 

import logging
import numpy as np
log = logging.getLogger(__name__)

class Computations:
    @classmethod
    def overall_mean(cls,values,mask):
        return np.average(values, weights=mask)
    
    @classmethod
    def overall_sigma(cls,values,mask):
        overall_mean_squared = np.average(values**2, weights=mask)
        return np.sqrt(overall_mean_squared - cls.overall_mean(values,mask)**2)

    @classmethod
    def mean_per_pixel(cls,values,mask=False):
        if not np.any(mask):
            return np.average(values,axis=0)
        else:
            w = mask.sum(axis=0)>0 # This is because one value is not masked I should put the axis value
            return np.average(values[:,w],weights=mask[:,w],axis=0)
    
    @classmethod
    def std_per_pixel(cls,values,mask):
        w = mask.sum(axis=0)>0
        mean = cls.mean_per_pixel(values,mask)
        mean_squared = cls.mean_per_pixel(values**2,mask)
        return np.sqrt(mean_squared-mean**2)

