function [] = plotMRmodel(model,magn,alpha);


fnms = fieldnames(model);
n_fld = size(fnms,1);
for jj = 1:n_fld
    
    clear cntnr.pts
    cntnr = getfield(model,fnms{jj});
    
    if strfind(fnms{jj},'clpCntr')
        continue
    end
    
    if strfind(fnms{jj},'clips')
        
        n_clips = size(cntnr,2);
        
        for ii = 1:n_clips
            plot3(cntnr(ii).pts(1,:)*magn,cntnr(ii).pts(2,:)*magn,cntnr(ii).pts(3,:)*magn,'Color',[0.6 0.6 0.6],'linewidth',1.3);
               text(cntnr(ii).pts(1,1)*magn,cntnr(ii).pts(2,1)*magn,cntnr(ii).pts(3,1)*magn,num2str(ii),'FontSize',20,'Color','r')

        end
        continue
    end
    
    
    % the function assumes that all structures in contours are organised in
    % struct with 3D points insisde the struct in field pts
    
%     if strfind(fnms{jj},'probability')
%         %plot3(cntnr.pts(1,:)*magn,cntnr.pts(2,:)*magn,cntnr.pts(3,:)*magn,'.g');
%         
%         continue
%     end
%     
    if strfind(fnms{jj},'optic_nerve')
        plot3(cntnr.pts(1,:)*magn,cntnr.pts(2,:)*magn,cntnr.pts(3,:)*magn,'Color',[0.5529    0.3608    0.9294],'linewidth',1);
        continue
        
    end
    
    if strfind(fnms{jj},'lens')
        plot3(cntnr.pts(1,:)*magn,cntnr.pts(2,:)*magn,cntnr.pts(3,:)*magn,'Color',[0.9294    0.1    0.1],'linewidth',1);
        continue
        
    end
    if strfind(fnms{jj},'ONH')
        plot3(cntnr.pts(1,:)*magn,cntnr.pts(2,:)*magn,cntnr.pts(3,:)*magn,'Color',[0.5922    0.4196    0.7882] ,'linewidth',1);
        continue
        
    end
%     
%     if strfind(fnms{jj},'sclera_drr')
% 
% %         n_points = length(cntnr.pts);
% %         
% %         for ii = 1:n_points
% %             plot3(cntnr.pts(1,ii),cntnr.pts(2,ii),cntnr.pts(3,ii),'+k','MarkerSize', 20)%,...
% %               %  'MarkerFaceColor',[cntnr.gr_lvls(ii) cntnr.gr_lvls(ii) cntnr.gr_lvls(ii)],'MarkerEdgeColor',[cntnr.gr_lvls(ii) cntnr.gr_lvls(ii) cntnr.gr_lvls(ii)]);
% %        
% %         end
% %         continue
%     end
%     
    if strfind(fnms{jj},'sclera')
        plot3(cntnr.pts(1,:),cntnr.pts(2,:),cntnr.pts(3,:),'ob','MarkerSize', 1);
        %                plot3(cntnr.pts(1,:)*magn,cntnr.pts(2,:)*magn,cntnr.pts(3,:)*magn,'Color','b','linewidth',0.4);
        %                plot3(cntnr.pts(1,:)*magn,cntnr.pts(2,:)*magn,cntnr.pts(3,:)*magn,'Color',[jj/n_fld 0 (n_fld-jj)/n_fld alpha],'linewidth',0.7);
        
        
    end
%     
%     if strfind(fnms{jj},'trgtBase')
%         plot3(cntnr.pts(1,:),cntnr.pts(2,:),cntnr.pts(3,:),'dk','MarkerSize', 3);
%         %                plot3(cntnr.pts(1,:)*magn,cntnr.pts(2,:)*magn,cntnr.pts(3,:)*magn,'Color','b','linewidth',0.4);
%         %                plot3(cntnr.pts(1,:)*magn,cntnr.pts(2,:)*magn,cntnr.pts(3,:)*magn,'Color',[jj/n_fld 0 (n_fld-jj)/n_fld alpha],'linewidth',0.7);
%         
%         
%     end
% 
%     
% %     
%     if strfind(fnms{jj},'macula')
%         % plot3(cntnr.pts(1,:)*magn,cntnr.pts(2,:)*magn,cntnr.pts(3,:)*magn,'Color','g','linewidth',1);
%         %plot3(cntnr.pts(1,:)*magn,cntnr.pts(2,:)*magn,cntnr.pts(3,:)*magn,'g','linewidth',1);
%         
%         %continue
%     end
%     

  if strfind(fnms{jj},'foveaEPTN')
        plot3(cntnr.pts(1,:)*magn,cntnr.pts(2,:)*magn,cntnr.pts(3,:)*magn,'ok','linewidth',1);
        %plot3(cntnr.pts(1,:)*magn,cntnr.pts(2,:)*magn,cntnr.pts(3,:)*magn,'g','linewidth',1);
        %continue
    end

  if strfind(fnms{jj},'maculaEPTN')
        plot3(cntnr.pts(1,:)*magn,cntnr.pts(2,:)*magn,cntnr.pts(3,:)*magn,'or','linewidth',1);
        %plot3(cntnr.pts(1,:)*magn,cntnr.pts(2,:)*magn,cntnr.pts(3,:)*magn,'g','linewidth',1);
        %continue
    end

    if strfind(fnms{jj},'macPrbbltyArea')
        plot3(cntnr.pts(1,:)*magn,cntnr.pts(2,:)*magn,cntnr.pts(3,:)*magn,'-og','linewidth',1);
        %plot3(cntnr.pts(1,:)*magn,cntnr.pts(2,:)*magn,cntnr.pts(3,:)*magn,'g','linewidth',1);
        %continue
    end
%     
    if not(isempty(strfind(fnms{jj},'tumor'))) || not(isempty(strfind(fnms{jj},'GTV')))
        %
        %           aa=makesurfaceMR(model.T1_tumor',0);
        %         trisurf(aa.surf.faces,aa.PP(strfind(fnms{jj},'tumor'):,1),aa.PP(:,2),aa.PP(:,3));
        plot3(cntnr.pts(1,:)*magn,cntnr.pts(2,:)*magn,cntnr.pts(3,:)*magn,'r','linewidth',1.3);
        continue
    end
   
    if strfind(fnms{jj},'tumEnl')
        %
        %           aa=makesurfaceMR(model.T1_tumor',0);
        %         trisurf(aa.surf.faces,aa.PP(:,1),aa.PP(:,2),aa.PP(:,3));
        plot3(cntnr.pts(1,:)*magn,cntnr.pts(2,:)*magn,cntnr.pts(3,:)*magn,'+k','linewidth',1.3);
        continue
    end
    
    if strfind(fnms{jj},'optcDsk')
        %
        %           aa=makesurfaceMR(model.T1_tumor',0);
        %         trisurf(aa.surf.faces,aa.PP(:,1),aa.PP(:,2),aa.PP(:,3));
        plot3(cntnr.pts(1,:)*magn,cntnr.pts(2,:)*magn,cntnr.pts(3,:)*magn,'-or','MarkerSize',10);
        continue
    end
%     
%     if strfind(fnms{jj},'target_thickness')
%         
%         continue
%     end
%     
    
    
    
    
end

end