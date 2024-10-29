function [rmodel] = readMRModelfromStrct(cont)

num_str = size(cont,2)

num_of_clips = 1;
for ii = 1:num_str
    if contains(cont(ii).ROIName,'clip') || contains(cont(ii).ROIName,'Clip')
        rmodel.clips(num_of_clips).pts = [cont(ii).Points'; ones(1,length(cont(ii).Points))];
        rmodel.clpCntr(num_of_clips,:) = mean(rmodel.clips(num_of_clips).pts(1:3,:)');
        num_of_clips = num_of_clips+1;
    else
        
        rmodel.(genvarname(cont(ii).ROIName)) = [cont(ii).Points'; ones(1,length(cont(ii).Points))];
        
    end
end


end