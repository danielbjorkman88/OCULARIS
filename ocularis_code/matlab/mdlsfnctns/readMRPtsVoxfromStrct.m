function [rmodel] = readMRPtsVoxfromStrct(cont)

num_str = size(cont,2)

num_of_clips = 1;
for ii = 1:num_str 
    if contains(cont(ii).ROIName,'clip') || contains(cont(ii).ROIName,'Clip')
        rmodel.clips(num_of_clips).pts = [cont(ii).Points'; ones(1,length(cont(ii).Points))];
        rmodel.clpCntr(num_of_clips,:) = mean(rmodel.clips(num_of_clips).pts(1:3,:)');
        num_of_clips = num_of_clips+1;
    else
        
        rmodel.(genvarname(cont(ii).ROIName)).pts = [cont(ii).Points'; ones(1,length(cont(ii).Points))];
        rmodel.(genvarname(cont(ii).ROIName)).vxp = [cont(ii).VoxPoints'; ones(1,length(cont(ii).VoxPoints))];
    
        
        n_sl = size(cont(ii).Segmentation,3);
        
        tot = [];
        for jj = 1:n_sl
           se = strel('disk',3);
           afterOpening = imopen(cont(ii).Segmentation(:,:,jj),se); 
           bw = edge(afterOpening);
           idxList = regionprops(bw,'PixelList');
           if isempty(idxList)
               continue
           else
               clear added
               added = [idxList(1).PixelList'; jj*ones(1,length(idxList(1).PixelList));ones(1,length(idxList(1).PixelList))];
               tot = cat(2,tot,added);
               
           end
            
        end
        
        rmodel.(genvarname(cont(ii).ROIName)).vxpBB = tot;
        rmodel.(genvarname(cont(ii).ROIName)).sgmVol = cont(ii).Segmentation;
    end
    
end








end