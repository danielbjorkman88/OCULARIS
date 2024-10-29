function [out] = deConTumourEP(points3D,plot)

out.surf.faces =[];
out.surf.vertices =points3D;
out.PP = points3D;
out.V = 0;

points3D(end+1,:) = points3D(1,:);

% indices to unique values in points3D
[c, ind] = unique(points3D, 'rows');
% duplicate indices
duplicate_ind = setdiff(1:size(points3D, 1), ind);

index = duplicate_ind(1:2:end);

[m,n] = size(index);

clear pp K V 
pp = points3D(1:index(2),:);

[~,V] = convhull(pp(:,1),pp(:,2),pp(:,3));
out.V = out.V +V;

out.surf.vertices = cat(1,out.surf.vertices,pp);

for tt = 1:index(2)-1
    if tt < index(1)
        K(tt,:) = [tt,tt+index(1),tt+1];
    else
        K(tt,:) = [tt,tt-index(1)+1,tt+1];
    end
end

out.surf.faces = cat(1,out.surf.faces,K);

if plot
close all
figure(1), hold on,
trisurf(K,pp(:,1),pp(:,2),pp(:,3))
end

for i = 2:n-1
    
    clear pp K V 
    pp = points3D(index(i-1):index(i+1),:);
    
    [~,V] = convhull(pp(:,1),pp(:,2),pp(:,3));
    out.V = out.V +V;
    
    ii=0;
    for tt = index(i-1):index(i+1)-1
        ii = ii + 1;
        if tt < index(i)
            K(ii,:) = [tt,tt+index(i)-index(i-1),tt+1];
        else
            K(ii,:) = [tt,tt-index(i)+1+index(i-1),tt+1];
        end
    end
    
    out.surf.faces = cat(1,out.surf.faces,K);
    
    if plot
        cla
        trisurf(out.surf.faces,points3D(:,1),points3D(:,2),points3D(:,3))
    end
end

clear pp K V 
pp = [points3D(index(i):index(i+1),:);
      points3D(1:index(1),:)];

[~,V] = convhull(pp(:,1),pp(:,2),pp(:,3));
out.V = out.V +V;

ii = 0;
for tt = index(i):index(i+1)+index(1)-2
    ii = ii + 1;
    if tt < index(i+1)
        K(ii,:) = [tt,tt-index(i)+1,tt+1];
    else
        K(ii,:) = [tt-index(i+1)+1,tt-index(i)+1+index(i-1),tt-index(i+1)+2];
    end
end
    



out.surf.faces = cat(1,out.surf.faces,K);


if plot
    cla
    trisurf(out.surf.faces,points3D(:,1),points3D(:,2),points3D(:,3))
%     plot3(points3D(:,1),points3D(:,2),points3D(:,3),'k','linewidth',10)
end

end