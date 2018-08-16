classdef information_analysis_functions < handle 
    methods (Static)
        %% histogram funcions
        
        function res = return_smoothed_hist_1D(rdat_,bw_,smw_)
            res = struct();
            
            rdat_nanpos = isnan(rdat_); 
            % create edges and centers for histogramming from non-nan data
            dat_nonnan = rdat_(~rdat_nanpos); 
            edge_dat = (min(dat_nonnan)-2*std(dat_nonnan)-bw_):bw_:(max(dat_nonnan)+2*std(dat_nonnan)+bw_);
            cent_dat = (edge_dat(1:end-1)+edge_dat(2:end))/2;
            
            % then replace nan data with largest center  
            dat_ = rdat_; 
            dat_(rdat_nanpos) = cent_dat(end);
            
            if ~ischar(smw_)  
                hist_dat = histcounts(dat_, edge_dat);
                sm_hist_dat = conv(hist_dat, smw_, 'same');
            else 
                if strcmpi(smw_, 'ksdensity') 
                    [sm_hist_dat,~] = ksdensity(dat_, cent_dat); 
                else
                    error('information_analysis_functions::return_smoothed_hist_1D: wrong ''smw_'' input');
                end
            end
            
            res.centers = cent_dat;
            res.edges = edge_dat;
            res.smoothed_hist = sm_hist_dat; 
            res.nonnan_data = dat_; 
            
            sm_hist_dat_normz = sm_hist_dat/sum(sm_hist_dat(:)); 
            res.entropy = information_analysis_functions.return_entropy(sm_hist_dat_normz,bw_); 
        end
        
        function res = return_smoothed_hist_2D_and_MI(dat1_,bw1_,dat2_,bw2_,sm_1D,sm_2D)
            res = struct();
            sm_h1 = information_analysis_functions.return_smoothed_hist_1D(dat1_,bw1_,sm_1D);
            sm_h2 = information_analysis_functions.return_smoothed_hist_1D(dat2_,bw2_,sm_1D);
            
            if ~ischar(sm_2D) 
                hist_joint = histcounts2(sm_h1.nonnan_data, sm_h2.nonnan_data, sm_h1.edges, sm_h2.edges);
                sm_hjoint = conv2(hist_joint, sm_2D, 'same');
            else 
                if strcmpi(sm_2D, 'ksdensity') 
                    [cent_1, cent_2] = meshgrid(sm_h1.centers, sm_h2.centers); 
                    cent_1 = cent_1(:); 
                    cent_2 = cent_2(:); 
                    cent_jnt = [cent_1, cent_2]; 
                    data_jnt = [transpose(sm_h1.nonnan_data), transpose(sm_h2.nonnan_data)]; 
                    [sm_hjoint,~] = ksdensity(data_jnt, cent_jnt); 
                    sm_hjoint = reshape(sm_hjoint, [length(sm_h2.centers),length(sm_h1.centers)]); 
                    sm_hjoint = transpose(sm_hjoint); 
                else
                    error('information_analysis_functions::return_smoothed_hist_2D_and_MI: wrong ''sm_2D'' input'); 
                end
            end
            
            sm_hjoint_normz = sm_hjoint/sum(sm_hjoint(:));             
            joint_entropy = information_analysis_functions.return_entropy(sm_hjoint_normz,[bw1_,bw2_]); 
            
            res.dat1 = sm_h1;
            res.dat2 = sm_h2;
            
            res.joint = struct('smoothed_hist', sm_hjoint, 'entropy', joint_entropy); 
            res.MI = sm_h1.entropy + sm_h2.entropy - joint_entropy;
        end
        
        %% entropy function
        function res = return_entropy(hist_,bw_)
            hist_ = hist_/sum(hist_(:)); % just in case not normz
            res = -hist_.*log2(hist_);
            res(~isfinite(res)) = 0;
            res = sum(res(:)) + 0*log2(prod(bw_));
        end
        
    end
    
end