function [rmodel,geom_str,dd_lens_eye] = readEPModelfromPln(PathName,FileName)


filename = [PathName FileName];

disp(filename)

Eye.Post_Pole{1} = 'Posterior Pole';
Eye.Post_Pole{2} = read_pln(filename, Eye.Post_Pole{1});
Eye.Macula{1} = 'Macula';
Eye.Macula{2} = read_pln(filename, Eye.Macula{1});
Eye.Disc{1} = 'Disc';
Eye.Disc{2} = read_pln(filename, Eye.Disc{1});
Eye.Limbus{1} = 'Limbus';
Eye.Limbus{2} = read_pln(filename, Eye.Limbus{1});
Eye.R_Clock{1} = 'Retina Clock'; % 3 times
Eye.R_Clock{2} = read_pln(filename, Eye.R_Clock{1});
Tumor.Surf{1} = 'Surface Of Tumor';
Tumor.Surf{2} = read_pln(filename, Tumor.Surf{1});
Tumor.Circ{1} = 'Circumference Of Tumor';
Tumor.Circ{2} = read_pln(filename, Tumor.Circ{1});
Lens.Surf{1} = 'Lens Surface';
Lens.Surf{2} = read_pln(filename, Lens.Surf{1});
os.Surf{1} = 'Ora Serrata';
os.Surf{2} = read_pln(filename, os.Surf{1});
Clip.C_File{1} = 'Clip File';
Clip.C_File{2} = read_pln(filename, Clip.C_File{1});

% Eye.Upper_lid{1} = 'Upper Lid';
% Eye.Upper_lid{2} = read_pln(filename, Eye.Upper_lid{1});

aoe.Surf{1} = 'Axis Of Eye';
aoe.Surf{2} = read_pln(filename, aoe.Surf{1});
%% cornea creation

r_g1 = abs(min(Eye.Post_Pole{2}(:,3)));
eye_length = r_g1/0.485;
r_g2 = eye_length*0.375;

z_g2 = eye_length - r_g1 -r_g2;

a_g1 = max([max(Eye.R_Clock{2}(:,1)) max(Eye.R_Clock{2}(:,2))]);
c_g1 = r_g1;

% anterior_boarder = Eye.Iris{2}(1,3);
anterior_boarder = Eye.Limbus{2}(1,3);

% cornea
[x_c,y_c,z_c] = ellipsoid(0,0,z_g2,r_g2,r_g2,r_g2);
mask_anterior = z_c >= anterior_boarder;
cornea = [x_c(mask_anterior) y_c(mask_anterior) z_c(mask_anterior);
    Eye.Limbus{2}(:,:)];
% main globe
[x_g,y_g,z_g] = ellipsoid(0,0,0,a_g1,a_g1,c_g1);
mask_main = z_g < anterior_boarder;
main_globe = [x_g(mask_main) y_g(mask_main) z_g(mask_main);
    Eye.Limbus{2}(:,:)];

whole_eye = [main_globe; cornea];

geom_str = cell(5,1);
geom_str{1} = [Tumor.Surf{2}(:,:) ones(size(Tumor.Surf{2}(:,1),1),1)]';
geom_str{2} = [whole_eye ones(size(whole_eye(:,1),1),1)]';
geom_str{3} = [cornea ones(size(cornea(:,1),1),1)]';
geom_str{4} = [Eye.Macula{2}(:,:) ones(size(Eye.Macula{2}(:,1),1),1)]';
geom_str{5} = [Eye.Disc{2}(:,:) ones(size(Eye.Disc{2}(:,1),1),1)]';
geom_str{6} = [Clip.C_File{2}(:,:) ones(size(Clip.C_File{2}(:,1),1),1)]';
geom_str{7} = [Lens.Surf{2}(:,:) ones(size(Lens.Surf{2}(:,1),1),1)]';
geom_str{8} = [Tumor.Circ{2}(:,:) ones(size(Tumor.Circ{2}(:,1),1),1)]';
geom_str{9} = [main_globe ones(size(main_globe(:,1),1),1)]';
% geom_str{10} = [Eye.Upper_lid{2}(:,:) ones(size(Eye.Upper_lid{2}(:,1),1),1)]';
axisofeye = [aoe.Surf{2}(:,:) ones(size(aoe.Surf{2}(:,1),1),1)]';

num_of_clips = size(geom_str{6},2)/20;

Acommon = intersect(whole_eye,cornea,'rows');
sclera= setxor(whole_eye,Acommon, 'rows');
[we_center,we_radius] = sphereFit(sclera);

pnts_lens = geom_str{7}(1:3,:)';
lens_center = fitEllipse3d(pnts_lens);
lens_center(4:end) = [];

clear idx
[mdMR,idx] = sort(pnts_lens,3);
z_ints = mean(mdMR(1:3,3)) 
[aa] = find(pnts_lens(:,3) == z_ints);

lens_apex = fitCircle3d(pnts_lens(aa,:));
lens_apex(4:end) = [];

pp_ = [geom_str{8}(1,:)',geom_str{8}(2,:)',geom_str{8}(3,:)'];

for jj = 1:size(pp_,1)
    dir = (pp_(jj,:)-we_center)./norm(pp_(jj,:)-we_center);
    pp(jj,:) = dir.*we_radius;
    
end
[az,elev,r] = cart2sph(pp(:,1),pp(:,2),pp(:,3));
s = referenceSphere('Unit Sphere','mm');
s.Radius = we_radius;
uu = areaint(rad2deg(elev),rad2deg(az),s);

rmodel.scleralR = we_radius;
rmodel.trgt_baseArea = uu;
rmodel.trgt_base = [pp';ones(1,length(pp))];
rmodel.lens_apex = lens_apex';
rmodel.eye_center = we_center';
dd_lens_eye = norm(lens_apex - we_center);

for jj = 1:num_of_clips
    rmodel.clips(jj).pts = geom_str{6}(:,1+(jj-1)*20:jj*20);
    
    rmodel.clpCntr(:,jj) = mean(rmodel.clips(jj).pts(1:3,:)');
end

rmodel.target = geom_str{1};
rmodel.eyeglobe = geom_str{2};
rmodel.cornea = geom_str{3};
rmodel.macula = geom_str{4};
rmodel.optdisk = geom_str{5};
rmodel.main_globe = geom_str{9};
rmodel.sclera = [sclera';ones(1,length(sclera))];
rmodel.axisofeye = axisofeye;
rmodel.cor = mean(axisofeye(1:3,:)')';
rmodel.lens = geom_str{7};
% rmodel.upper_lid = geom_str{10};
% 3D disc
vec_1 = [geom_str{5}(1,1),geom_str{5}(2,1),geom_str{5}(3,1)] -[geom_str{5}(1,end),geom_str{5}(2,end),geom_str{5}(3,end)];
vec_2 = [geom_str{5}(1,2),geom_str{5}(2,2),geom_str{5}(3,2)] -[geom_str{5}(1,end),geom_str{5}(2,end),geom_str{5}(3,end)];
orto_norm = cross(vec_1,vec_2)./sqrt(dot(cross(vec_1,vec_2),cross(vec_1,vec_2)));
thickness_disc = 0.05;
disc_3D = [Eye.Disc{2}(:,:);
    Eye.Disc{2}(:,:) + repmat(orto_norm.*thickness_disc,size(Eye.Disc{2}(:,:),1),1)];
geom_str{5} = [disc_3D ones(size(disc_3D(:,1),1),1)]';


abs_volume = zeros(size(geom_str));
for i = 1:length(geom_str)
    structure_i = geom_str{i}';
    DT = delaunayTriangulation(structure_i(:,1),structure_i(:,2),structure_i(:,3));
    [FBtri,FBpoints] = freeBoundary(DT); % only surface of object, no tetraedons
    TR_Polygon = triangulation(FBtri,FBpoints);
    abs_volume(i) = Vol_Polygon_pyramids(TR_Polygon);
end

display(['volume of the eye: ' num2str(abs_volume(2)) ' mm^3'])