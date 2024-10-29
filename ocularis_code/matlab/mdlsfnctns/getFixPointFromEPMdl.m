function [pA,pt3d] = getFixPointFromEPMdl(mdl)

los = createLine3d(mdl.eye_center',mdl.lens_apex');
planeFix = createPlane([0, 0, 132.5],[0,0,-1]);
pt3d = intersectLinePlane(los,planeFix);
pA = getFixPntPA(pt3d);

end