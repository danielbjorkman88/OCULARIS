function [] = plotEPeyeGlobe(eyeGlobePointR,magn,alpha)

[~,idx]=unique(eyeGlobePointR,'rows','first');
uu=eyeGlobePointR(sort(idx),:);
for jj =1:15:size(uu,1)-170
    uu(jj,:) = [];
    rng = jj:jj+15;
    p(1) = plot3(uu(rng,1)*magn,uu(rng,2)*magn,uu(rng,3)*magn,'Color',[0 .13 .78 alpha],'linewidth',2);

end

finRng1 = rng(end)+1:rng(end)+60;
 p(1) = plot3(uu(finRng1,1)*magn,uu(finRng1,2)*magn,uu(finRng1,3)*magn,'Color',[0 .13 .78 alpha],'linewidth',2);
finRng2 = finRng1(end)+1:length(uu);
 p(1) = plot3(uu(finRng2,1)*magn,uu(finRng2,2)*magn,uu(finRng2,3)*magn,'Color',[0 .13 .78 alpha],'linewidth',2);

 


end