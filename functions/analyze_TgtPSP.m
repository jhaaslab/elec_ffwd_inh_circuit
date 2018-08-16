function result = analyze_TgtPSP(tv, vm)
psp = vm - vm(1); 
if isempty(psp(psp>0)) 
    peak_epsp = 0;
    auc_epsp = 0; 
    dt_epsp = 0;
else
    peak_epsp = max(psp); 
    t_peak_epsp = tv(psp == peak_epsp); 
    if min(psp) < 0 % ipsp exists
        upperbound_2findtepsp = tv(psp == min(psp)); 
    else % pseudo-tpeak-ipsp to find int. window before this time point 
        % actually the time point of 5% peak_epsp after tpeak-epsp 
        upperbound_2findtepsp = min(tv(tv >= t_peak_epsp & psp <= 0.05*peak_epsp)); 
    end 
    loc_epsp = psp > 0 & tv < upperbound_2findtepsp; % avoiding ADP 
    epsp_rect = zeros(size(tv)); 
    epsp_rect(loc_epsp) = 1; 
    
    dt_epsp = trapz(tv, epsp_rect);
    
    epsp_vec = psp.*epsp_rect; 
    auc_epsp = trapz(tv, epsp_vec); 
    
        
end
result = struct('peak_epsp', peak_epsp, 'dt_epsp', dt_epsp, 'auc_epsp', auc_epsp);  
end
