function [MRmodel] = readMrModelFromDCM(fldr_path,img_path)

if nargin < 2
    imagedir = [fldr_path 'imgs/'];
else
    imagedir = img_path;
end


% output MRI model in mm space 
% input path to folder conainting exclusively dcm structures and images in
% subfolder imgs

dataPattern = fullfile(fldr_path,'*.dcm');
data_list = dir(dataPattern);
[foo,cont] = dicomrt2matlab([fldr_path data_list(1).name],imagedir);
[MRmodel] = readMRModelfromStrct(cont);



end
