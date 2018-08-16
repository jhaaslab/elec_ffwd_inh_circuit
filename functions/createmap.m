function cmap = createmap(color_options, ncolors, bound_tol, bound_darklvl)

if ischar(color_options) 
    switch lower(color_options) 
        case 'bwr'
            csamp_mat = [ 0,0,1; 
                          1,1,1;
                          1,0,0]; 
        case 'bwg' 
            csamp_mat = [ 0,0,1; 
                          1,1,1;
                          0,0,1]; 
        case 'pwo'
            csamp_mat = [ [102,0,255]/255;
                          [1,1,1]*0.98; 
                          [255,102,0]/255];
        otherwise 
            error(['createmap::''color_options'' can only be a string like '...
                '''bwr'', ''bwg'' or ''pwo'', or a cell of colors, not %s'], color_options);
    end
else
    if iscell(color_options) 
        csamp_mat = vertcat(color_options{:}); 
    else
        error(['createmap::''color_options'' can only be a string '...
            ' or a cell of colors']);
    end
end
nsamp = size(csamp_mat,1); 
samp_vec = linspace(0,1,nsamp);

if strcmpi(ncolors, 'auto')
    ncolors = size(get(gcf,'colormap'),1); 
end
interp_vec = linspace(0,1,ncolors); 
cinterp_mat = zeros(ncolors,3);

rect_ = ones(size(interp_vec))*0.8; 
rect_(interp_vec > bound_tol & interp_vec < (1-bound_tol)) = 1; 
hannwind = hanning(ceil(length(rect_)/10));
hannwind = hannwind/sum(hannwind(:)); 
rect_ = conv(rect_, hannwind, 'same');  
rect_ = (1-bound_darklvl)*(rect_-min(rect_))/(max(rect_)-min(rect_)) + bound_darklvl; 

for i = 1:3 
    cinterp_mat(:,i) = rect_.*interp1(samp_vec, csamp_mat(:,i), interp_vec, 'pchip');
end 

cinterp_mat = cinterp_mat/max(cinterp_mat(:)); 

cmap = cinterp_mat; 


end