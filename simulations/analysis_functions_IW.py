from __future__ import division
import numpy as np 

def calc_IW(tv, vm): 
    psp = vm - vm[0] # compared to baseline - considered to be v(0) 

    if np.size(psp[psp>0]) == 0:
        peak_epsp = 0 
        auc_epsp = 0 
        dt_epsp = 0
    else: 
        peak_epsp = max(psp)
        t_epsp_max = tv[psp == max(psp)] 
        if min(psp) < 0:
            t_ipsp_max = tv[psp == min(psp)][0] 
        else: 
            t_ipsp_max = min(tv[np.logical_and(tv>=t_epsp_max, psp <= 0.05*peak_epsp)]) # 5% of epsp_peak if does not have ipsp 
        loc_epsp = np.logical_and(psp > 0, tv<t_ipsp_max) # to watch out for situations like ADP after being hyperpolarized
        
        dt_rect = np.zeros(np.size(tv))
        dt_rect[loc_epsp] = 1 
        dt_epsp = np.trapz(dt_rect ,tv) #max(t_epsp) - min(t_epsp) only if known 2nd epsp is not big enough 
        
        auc_rect = np.zeros(np.size(psp)) 
        auc_rect[loc_epsp] = psp[loc_epsp]
        auc_epsp = np.trapz(auc_rect, tv)

    return {'peak_epsp': peak_epsp, 'dt_epsp': dt_epsp, 'auc_epsp': auc_epsp}
