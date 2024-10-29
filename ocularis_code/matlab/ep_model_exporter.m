%% Reading and geometrical sorting of EyePlan Dose Data
clc,clear all,close all
addpath('mdlsfnctns/')
addpath('geom3d/')



base_path = 'P:\OPTIS2\OptisTPS\patient_database';

base_path = 'Z:\PatientData';

pln_file = '12345PATIENTNAME_Pln_4.pln';
xpc_fn = '12345PATIENTNAME_4.xpc';



patient_number = pln_file(1:5);


read_path = append(base_path , '\P' , patient_number , '/');


write_path = base_path + "\P" + patient_number + "\" + 'model/';
write_path_registered = base_path + "\P" + patient_number + "\" + 'model_clips/';


transform_matrix_filename = append(patient_number , '_' , pln_file(end-4) , '_clips_transformation.csv');


if ~exist(write_path, 'dir')
    % If it doesn't exist, create it
    mkdir(write_path);
    disp(['Directory created: ', write_path]);
else
    disp(['Directory already exists: ', write_path]);
end




%% Read Eye Plan Model
modEP = readEPModelfromPln(read_path,pln_file); %modEP eye structures in gaze centered position


% Writes the content of modEP, ie the structures in gaze centered
%reference frame
fn = fieldnames(modEP);
for k=1:numel(fn)
    if( isnumeric(modEP.(fn{k})) )
        disp(fn{k})
        value = getfield(modEP,fn{k});
        out_filename = pln_file(1:5) + "_" + pln_file(length(pln_file) -4) +"_clips" + "_" + fn{k} + ".csv";
        dlmwrite(fullfile(write_path, out_filename), value, 'precision', 12);
    end
end

% return


if ~exist(write_path_registered, 'dir')
    % If it doesn't exist, create it
    mkdir(write_path_registered);
    disp(['Directory created: ', write_path_registered]);
else
    disp(['Directory already exists: ', write_path_registered]);
end

%% Registration using clips
[axial,lat,num_of_clips] = gClipsFromFile([read_path xpc_fn]);
clips3D = triangClips(axial,lat,num_of_clips);
[r,t] = globale_new(modEP.clpCntr',clips3D);
T = costrrmat(r); % T = Transformation matrix
modEPf = trnsfrmEPmodel(modEP,T); %modEPf eye structures in treatment position



% Writes the transformation matrix
dlmwrite(fullfile(write_path_registered, transform_matrix_filename) , T, 'precision', 12);
 

% Writes the content of modEPf, ie the structures in treatment orientation 
fn = fieldnames(modEPf);
for k=1:numel(fn)
    if( isnumeric(modEPf.(fn{k})) )
        disp(fn{k})
        value = getfield(modEPf,fn{k});
        out_filename = pln_file(1:5) + "_" + pln_file(length(pln_file) -4) +"_clips" + "_" + fn{k} + ".csv";
        dlmwrite(fullfile(write_path_registered, out_filename), value, 'precision', 12);
    end
end



