function [modout] = harmoniseMRItoEPmodel(modelin)

fnms = fieldnames(modelin);
disp(fnms)
n_fld = size(fnms,1);
for jj = 1:n_fld
    clear cntnr
    cntr = getfield(modelin,fnms{jj});
    
    strct_name = lower(fnms{jj});
    
    if contains(strct_name,'cornea')
        modelout.cornea = cntr;
    elseif contains(strct_name,'lens')
        modelout.lens = cntr;
    elseif contains(strct_name,'eyeglobe') || contains(strct_name,'sklera') || contains(strct_name,'sclera')
        modelout.eyeglobe = cntr;
        %     elseif contains(strct_name,'sklera') || contains(strct_name,'sclera')
        %         modelout.sclera = cntr;
    elseif contains(strct_name,'gtv') || contains(strct_name,'tumor')
        modelout.target = cntr;
    elseif contains(strct_name,'clips')
        modelout.clips = cntr;
    elseif contains(strct_name,'clpcntr')
        modelout.clpCntr = cntr';
    elseif contains(strct_name,'Eyelid') || contains(strct_name,'eyelid')
        modelout.eyelid = cntr;
    end
%     disp(strct_name)
%     disp(size(cntr,1));
%     disp(size(cntr,2));
end

[eye_C,r] = sphereFit(modelout.eyeglobe(1:3,:)');
disp(fieldnames(modelout))

%% model Clips as disk
for jj = 1:size(modelout.clips,2)
    clip_centroid =  modelout.clpCntr(:,jj);
    modelout.clips(jj).pts = [(circlePlane3D(clip_centroid,(clip_centroid-eye_C')./norm(clip_centroid-eye_C'),1.25,0.5))',clip_centroid];
end
    

%% extra stuff for Gaze Direction


[center, radii, evecs, pars ] = ellipsoid_fit_new( modelout.lens(1:3,:)' );
lensCenter = center';
line = createLine3d(eye_C,lensCenter);
sphere = [lensCenter,min(radii)];
pts = intersectLineSphere(line,sphere);

dd = 0;
for jj = 1:size(pts,1)
    
    dd_ = norm(pts(jj,:) - eye_C);
    if dd_ > dd
        modelout.lens_apex = pts(jj,:)';
        dd = dd_;
    end
end

modelout.cor = eye_C';
sphere = [eye_C,r] ;
pts = intersectLineSphere(line,sphere);

dd = (modelout.lens_apex' - eye_C)./norm(modelout.lens_apex' - eye_C);
for jj = 1:size(pts,1)
    dd_ = (pts(jj,:) - eye_C)./norm(pts(jj,:) - eye_C);
    if sign(dd_) == sign(-dd)
        opticdisk_centroid = pts(jj,:)';
    end
end

modelout.optdisk = [(circlePlane3D(opticdisk_centroid,(opticdisk_centroid-eye_C')./norm(opticdisk_centroid-eye_C'),0.75,0.5))',opticdisk_centroid];
%% copy macula at same place for now 
modelout.macula = [(circlePlane3D(opticdisk_centroid,(opticdisk_centroid-eye_C')./norm(opticdisk_centroid-eye_C'),0.75,0.5))',opticdisk_centroid];

modelout.axisofeye = [modelout.lens_apex,lensCenter',modelout.cor,opticdisk_centroid];


%% Center model
T_t = eye(4,4);
T_t(1:3,4) = -modelout.cor;
modelt = trnsfrmEPmodel(modelout,T_t);

%% Rotate
dd_axis = modelt.lens_apex' - modelt.cor';
vc = dd_axis./norm(dd_axis);
R2 = getRotMatrixfrom2Vec(vc,[0 0 1]);
T_r = eye(4,4);T_r(1:3,1:3) = R2;
modout = trnsfrmEPmodel(modelt,T_r);


%% Fitting ellipsoid for Sclera reconstruction
xyz = modout.eyeglobe';
A = xyz;
A(xyz(:,3)<0,:) = [];
B = modout.cornea';

k = convhull(B(:,1),B(:,2));
for jj = 1:size(A,1)
  bcenter = A(jj,1:2);
  in = inpolygon( bcenter(2) , bcenter(1) , B(k,1) , B(k,2) );
  dd(jj) = in;
end

A(logical(dd),:) = [];
A = [A;xyz(xyz(:,3)<0,:)];

[ center, radii, evecs, v, chi2 ] = ellipsoid_fit_new( A, '' );
 %draw fit
mind = min( [ A ] );
maxd = max( [ A ] );
nsteps = 50;
step = ( maxd - mind ) / nsteps;
%[ x, y, z ] = meshgrid( linspace( mind(1) - step(1), maxd(1) + step(1), nsteps ), linspace( mind(2) - step(2), maxd(2) + step(2), nsteps ), linspace( mind(3) - step(3), maxd(3) + step(3), nsteps ) );
steps = 0.3;
[x,y,z] = meshgrid(mind(1):steps:maxd(1),mind(2):steps:maxd(2),mind(3):steps:maxd(3));
Ellipsoid = v(1) *x.*x +   v(2) * y.*y + v(3) * z.*z + ...
2*v(4) *x.*y + 2*v(5)*x.*z + 2*v(6) * y.*z + ...
2*v(7) *x    + 2*v(8)*y    + 2*v(9) * z;
[faces,verts] =  isosurface( x, y, z, Ellipsoid, -v(10) );

verts(verts(:,3)<0,:) = [];

for jj = 1:size(verts,1)
  bcenter = verts(jj,1:2);
  in = inpolygon( bcenter(2) , bcenter(1) , B(k,1) , B(k,2) );
  dd(jj) = in;
end

verts(logical(dd),:) = [];

modout.sclera_resampled = verts';

%% Visualise
% figure,hold on,axis equal,
% plotEPmodel(modout,-1,1);
% p = patch( isosurface( x, y, z, Ellipsoid, -v(10) ) );


%% Old stuff (simple resamlipng with Z>0
% xyz = modout.eyeglobe;
% xyz(:,xyz(3,:)<0) = [];x = xyz(1,:);y = xyz(2,:);z = xyz(3,:);
% xlin = min(x):0.1:max(x);
% ylin = min(y):0.1:max(y);
% % xlin = linspace(min(x), max(x), 100);
% % ylin = linspace(min(y), max(y), 100);
% [X,Y] = meshgrid(xlin, ylin);
% Z = griddata(x,y,z,X,Y,'natural');
% x_ = reshape(X,1,[]);
% y_ = reshape(Y,1,[]);
% z_ = reshape(Z,1,[]);
% ii = isnan(z_);
% x_(ii) = [];
% y_(ii) = [];
% z_(ii) = [];
% 
% modout.eyeglobe_resampled = [x_;y_;z_];



end