function ePModTransf = trnsfrmEPmodel(model,M)

fnms = fieldnames(model);

for jj = 1:size(fnms)
    
    clear urca 
    urca = getfield(model,fnms{jj})';

    %strct_name = lower(fnms{jj});
    %disp(strct_name);
    
    
    if isstruct(urca)
       % For clips

       % it means that it is either the clips or another nested struct..
       % keep in mind that in ordr to work the nested variable must be
       % named pts
       int_fnms = size(urca,1);
       for tt = 1:int_fnms
           m.(genvarname(fnms{jj}))(tt).pts = mapping(urca(tt).pts',M)';
           
       end
    else
        ss = length(urca);
        %disp(size(urca,1));
        %disp(size(urca,2));
        if ss >1 
            m.(genvarname(fnms{jj})) = mapping(urca,M)';
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



