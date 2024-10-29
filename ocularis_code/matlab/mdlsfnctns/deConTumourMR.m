function [out] = deConTumourMR(points3D,plot)
out.surf.faces =[];
out.surf.vertices =[];
out.PP = points3D;
out.V = 0;

pl = fitPlane(points3D(1:3,:));

xx = pl(4:6); yy=pl(7:9);
zz = cross(xx,yy);


ptsTrs = transformPoint3d(points3D,[xx',yy',zz']');
[ptsOrd,idx] = sortrows(ptsTrs,3);%,'descend');

% indices to unique values in column 3
[index_] = find(abs(diff(ptsOrd(:,3)))>0.1);
aa = find(index_> size(ptsOrd,1)-4);
index_(aa) = [];
index = [index_ ;length(ptsOrd)];


[m,n] = size(index);

iiDx = idx(1:index(1));
pp3D(1:index(1),:) = points3D(sort(iiDx),:);
for jj = 1:m-1
    clear iiDx
    iiDx = idx(index(jj)+1:index(jj+1));
    pp3D(index(jj)+1:index(jj+1),:) = points3D(sort(iiDx),:);
end


clear pp K V 

for tt = 1:length(index)
    if index(tt) > 3
        st_idx = tt;
        break;
    end
end

if st_idx == 1
    st_idx = 2;
end

pp = pp3D(1:index(st_idx),:);

[K,V] = convhull(pp(:,1),pp(:,2),pp(:,3));
out.V = out.V +V;

% ii = 0;
% for tt = 1:index(1)-1
%     ii = ii +1;
%     K(ii,:) =  [1 tt+1 tt+2];
% end
% ii = ii+1;
% K(ii,:) =  [1 tt+1 2];

for jj = 1:length(K)
   xx = K(jj,:);
%    if all(xx >= 1 & xx <= index(1))
%        K(jj,:) = [NaN NaN NaN];
%    end
   if all(xx > index(1) & xx <= index(2))
       K(jj,:) = [NaN NaN NaN];
   end
end
K(~any(~isnan(K), 2),:)=[];
out.surf.faces = cat(1,out.surf.faces,K);

if plot
%     close all
 %   figure, hold on, axis equal
    plot3(points3D(:,1),points3D(:,2),points3D(:,3),'+k')

    trisurf(K,pp(:,1),pp(:,2),pp(:,3))
end

for i = st_idx:m-2
    
    clear pp K V 
    pp = pp3D(index(i-1)+1:index(i+1),:);
    
    [K,V] = convhull(pp(:,1),pp(:,2),pp(:,3));
    out.V = out.V +V;
    
    K = K + index(i-1);
    
    for jj = 1:length(K)
        xx = K(jj,:);
        if all(xx > index(i-1) & xx <= index(i))
            K(jj,:) = [NaN NaN NaN];
        end
        if all(xx > index(i) & xx <= index(i+1))
            K(jj,:) = [NaN NaN NaN];
        end
    end
    K(~any(~isnan(K), 2),:)=[];
    
    out.surf.faces = cat(1,out.surf.faces,K);
    
    if plot
        cla
        plot3(points3D(:,1),points3D(:,2),points3D(:,3),'+k')
        trisurf(K,pp3D(:,1),pp3D(:,2),pp3D(:,3))
    end
end


pp = pp3D(index(i)+1:end,:);

[K,V] = convhull(pp(:,1),pp(:,2),pp(:,3));
out.V = out.V +V;

K = K + index(i);

for jj = 1:length(K)
    xx = K(jj,:);
    if all(xx > index(i) & xx <= index(i+1))
        K(jj,:) = [NaN NaN NaN];
    end
%     if all(xx > index(i) & xx <= index(i+1))
%         K(jj,:) = [NaN NaN NaN];
%     end
end
K(~any(~isnan(K), 2),:)=[];

out.surf.faces = cat(1,out.surf.faces,K);

out.surf.vertices = pp3D;
out.surf.verticer(:,4) = 1;
if plot
    cla
    plot3(points3D(:,1),points3D(:,2),points3D(:,3),'k','linewidth',10)
    trisurf(out.surf.faces,...
        out.surf.vertices(:,1),out.surf.vertices(:,2),out.surf.vertices(:,3))
end



end