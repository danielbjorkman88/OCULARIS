
function [] = plotEPmodel(rmodel,magn,alpha)

  plotEPeyeGlobe(rmodel.eyeglobe',magn,alpha);
% % 
p(2) = plot3(rmodel.target(1,:)*magn,rmodel.target(2,:)*magn,rmodel.target(3,:)*magn,'Color',[1 0.25 0 alpha],'linewidth',3);
 p(4) = plot3(rmodel.optdisk(1,:)*magn,rmodel.optdisk(2,:)*magn,rmodel.optdisk(3,:)*magn,'Color',[.13 .78 .78 alpha],'linewidth',3);
 p(4) = plot3(rmodel.macula(1,:)*magn,rmodel.macula(2,:)*magn,rmodel.macula(3,:)*magn,'Color',[0 .78  0.13 alpha],'linewidth',3);
 p(6) = plot3(rmodel.lens(1,:)*magn,rmodel.lens(2,:)*magn,rmodel.lens(3,:)*magn,'Color',[.78 0 .13 alpha],'linewidth',2);
% plot3c(rmodel.trgt_base(1,:)*magn,rmodel.trgt_base(2,:)*magn,rmodel.trgt_base(3,:)*magn,'o','Color',[1 0 1 alpha],'linewidth',2);



for jj = 1:size(rmodel.clips,2)
    plot3c(rmodel.clips(jj).pts(1,2:end)*magn,rmodel.clips(jj).pts(2,2:end)*magn,rmodel.clips(jj).pts(3,2:end)*magn,...
        'Color',[.0 0 .0 alpha],'linewidth',2);
     text(rmodel.clpCntr(1,jj)+2,rmodel.clpCntr(2,jj)+2,rmodel.clpCntr(3,jj)+2,num2str(jj),'FontSize',20)

    patch(rmodel.clips(jj).pts(1,:)*magn,rmodel.clips(jj).pts(2,:)*magn,rmodel.clips(jj).pts(3,:)*magn,[0.8 0.8 0.8]);
end

end

function plot3c(x,y,z,varargin)  
    x = [x(:) ; x(1)];   
    y = [y(:) ; y(1)];
    z = [z(:) ; z(1)];
    plot3(x,y,z,varargin{:})  
end