function ePModTransf = rotEPmodelAroundPoint(model,M,cntr)

fnms = fieldnames(model);

T_ = eye(4,4);
T_(1:3,4) = -cntr;

T = eye(4,4);
T(1:3,4) = cntr;

for jj = 1:size(fnms)
    
    clear urca 
    urca = getfield(model,fnms{jj})';
    
    if isstruct(urca)
        
        
       % it means that it is either the clips or another nested struct..
       % keep in mind that in ordr to work the nested variable must be
       % named pts
       int_fnms = size(urca,1);
       for tt = 1:int_fnms
           
           mT_= mapping(urca(tt).pts',T_);
           mRT_ = mapping(mT_,M);
           m.(genvarname(fnms{jj}))(tt).pts = mapping(mRT_,T)';
           
       end
    else
        ss = length(urca);
        if ss >1 
           mT_= mapping(urca,T_);
           mRT_ = mapping(mT_,M);
           m.(genvarname(fnms{jj})) = mapping(mRT_,T)';
        else
            m.(genvarname(fnms{jj})) = urca;
        end
            
    end
        
    
end

% m.eyeGlobePointR = mapping(model.eyeGlobePointR,M);
% m.tumorPointR = mapping(model.tumorPointR,M);
% m.macula.pointR = mapping(model.macula.pointR,M);
% m.op_discPointR = mapping(model.op_discPointR,M);
% m.anteriorPointR = mapping(model.anteriorPointR,M);            
% m.clipsCentr = mapping(model.clipsCentr,M);
% m.lensPointR = mapping(model.lensPointR,M);
% 
% m.fitTipOfTheEye = mapping(model.fitTipOfTheEye,M);
% m.corfitCenter = mapping(model.corfitCenter,M);
% m.sclfitCenter = mapping(model.sclfitCenter,M);
% m.lenfitCenter = mapping(model.lenfitCenter,M);
% 
% m.lenfitRadius = model.lenfitRadius;
% m.sclfitRadius = model.sclfitRadius;

ePModTransf = m;
end



