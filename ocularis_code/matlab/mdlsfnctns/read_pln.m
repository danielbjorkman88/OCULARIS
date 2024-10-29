function pattern_structure = read_pln(filename, pattern)
% structures in the .pln file are searched by name (Apex of Tumor)
% the structure of the file itself is as followed:
% 1) Befor each string are 2 Bytes in little endian order, which give the
% length of the String.
% 2) After the searched string 2 more strings follow (also wicth length),
% Afterwards are 2 'dead' Bytes, no content
% 3) in the next 2 Bytes (little endian) is the length of the whole
% structure given
% 4) flot32, float32, float32,  16bit short, 
% length of string (little endian), string
% 5) read out in a for loop

%filename = '23161VATELLIS_Pln_4.pln';
fileID = fopen(filename);
A = fread(fileID);
fclose(fileID);
%%%%%%%%% Tumor %%%%
% pattern = 'Circumference Of Tumor';
% pattern = 'Apex Of Tumor';
% pattern = 'Surface Of Tumor'; % 4 times!!! different parts/tumor? 
% stored tumor?!
%%%%%%%%%% Clips %%%%
% pattern = 'Clip File';
%%%%%%%%% Eye %%%%
% pattern = 'Sclera';
% pattern = 'Retina';
% pattern = 'Axis Of Eye';
% pattern = 'Retina Clock'; % 3 times
% pattern = 'Iris Clock';
% pattern = 'Posterior Pole';
% pattern = 'Macula';
% pattern = 'Disc';
% pattern = 'Optic Nerve';
% pattern = 'Lens Periphery';
% pattern = 'Lens Surface';
% pattern = 'Iris';
% pattern = 'Limbus';
% pattern = 'Equator';
% pattern = 'Ora Serrata';
% pattern = 'Skin Plane';

c_string = char(A)';
k = strfind(c_string,pattern);
if strcmp(pattern,'Retina') || strcmp(pattern,'Sclera')
    k = k(1);
end
pattern_structure = [];
if isempty(k)
    display('String not found')
else
    for j = 1:length(k)
        %%%%%%%% Search with fread %%%%%%
        fileID = fopen(filename);
        skip1 = fread(fileID,k(j)-3);
        str1 = fread(fileID,1,'*uint16','l');
        skip2 = char(fread(fileID,str1)');
        str2 = fread(fileID,1,'*uint16','l');
        skip3 = char(fread(fileID,str2)');
        str3 = fread(fileID,1,'*uint16','l');
        skip4 = char(fread(fileID,str3)');
        leer = fread(fileID,1,'*uint16','l');
        structure_length = fread(fileID,1,'*uint16','l');

        point_store = zeros(structure_length,3);
        short_store = cell(structure_length,1);
        str_store = cell(structure_length,1);

        for i =1:structure_length
            point_store(i,:) = fread(fileID,3,'float32');
            short_store{i} = fread(fileID,1,'uint16');
            str4 =fread(fileID,1,'*uint16','l');
            str_store{i} = char(fread(fileID,str4)');
        end
        fclose(fileID);
%         display(point_store')
%         display(short_store')
%         display(str_store')
        pattern_structure = [pattern_structure; point_store];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%% Search in byte matrix %%%%%%%%%
% pln = dec2bin((A));
% buf_apex = k;
% str1 = bin2dec([pln(k-1,:) pln(k-2,:)]); % length of apex string
% buf_apex = buf_apex + str1 + 2; 
% str2 = bin2dec([pln(buf_apex-1,:) pln(buf_apex-2,:)]); % length of Tumor
% buf_apex = buf_apex + str2 + 2;
% str3 = bin2dec([pln(buf_apex-1,:) pln(buf_apex-2,:)]); % length of E
% buf_apex = buf_apex + str3 + 2;
% str4 = bin2dec([pln(buf_apex-1,:) pln(buf_apex-2,:)]) % length of E
% buf_apex = buf_apex + str4 + 2;
% apex_l = bin2dec([pln(buf_apex-1,:) pln(buf_apex-2,:)]); % length of E
% buf_apex = buf_apex + 2;
% 
% Apex = zeros(apex_l,3);
% short_store = cell(apex_l,1);
% srt_store = cell(apex_l,1);
% for i = 1:apex_l
%     [pln(buf_apex,:):pln(buf_apex+4,:)]
% 
% end
% bin_float = pln(buf_apex:buf_apex+3,:);
% bin_float = bin_float(:)';
% NumAsFloat32 = typecast(uint32(bin2dec(bin_float)),'float')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%