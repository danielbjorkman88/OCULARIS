function clips3D = triangClips(axial,lat,num_of_clips)


ddSourceAx = -1244;
ddFilmAx = 617;

ddSourceLat = 1215;
ddFilmLat = -609;

sourceLat = [ddSourceLat 0 0];
filmLat = [ddFilmLat 0 0];

sourceAx = [0 0 ddSourceAx];
filmAx = [0 0 ddFilmAx];


for jj = 1:num_of_clips
    clear aa dd 
    aa = fitzgibboellipse(axial(jj).xy(:,1:2));
    dd = convert_coeff2ellips(aa);
    axial_c(jj,:) = [dd(1:2) ddFilmAx];
   
%     axial(jj).xy(:,3) = ddFilmAx;
%     axial_c(jj,:) = mean(axial(jj).xy);
    
    clear aa dd 
    aa = fitzgibboellipse(lat(jj).xy(:,1:2));
    dd = convert_coeff2ellips(aa);
    lat_c(jj,:) = [ddFilmLat dd(1:2)];
end

for jj = 1:num_of_clips
    plane(jj,:) = createPlane(lat_c(jj,:),sourceLat,filmLat);
    lineAx(jj,:) = createLine3d(axial_c(jj,:),sourceAx);
    lineLat(jj,:) = createLine3d(lat_c(jj,:),sourceLat);
    
    clips3D(jj,:) = intersectPlaneLine(plane(jj,:), lineAx(jj,:));
end
end