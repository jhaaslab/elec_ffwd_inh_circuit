function [ptsh, centers] = return_histogram(spk_times, n_trials, nbins, edge_lim, smooth_win)
edges = linspace(edge_lim(1), edge_lim(2), nbins+1); 
centers = (edges(1:end-1)+edges(2:end))/2; 
counts = histcounts(spk_times, edges);
counts = counts/n_trials;
sm_wind = hanning(smooth_win); % gausswin(smooth_win,smooth_win/10);
ptsh = conv(counts, sm_wind, 'same')/sum(sm_wind);
end