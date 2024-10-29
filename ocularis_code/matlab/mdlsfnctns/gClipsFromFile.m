function [axial,lat,num_of_clips] = gClipsFromFile(filename)


fid = fopen(filename);

count_row =0;
for jj =1:10000
tline = fgets(fid);
while ischar(tline)
   count_row = count_row+1;
   break
end
end

fid = fopen(filename);
for int = 1:18
    tline = fgets(fid);
    if int == 17
        lineTxt = strsplit(tline);
        num_of_clips = str2num(lineTxt{1,size(lineTxt,2)-1});
    end
    if int == 18
        lineTxt = strsplit(tline);
        num_of_points = str2num(lineTxt{1,size(lineTxt,2)-1});
    end
end
fclose(fid);

axial = [];
lateral = []


hdr = 4;
fid = fopen(filename);
clips_coord = 0;
count_clips = 1;
count_points = 1;

for int = 1:count_row
    tline = fgets(fid);
    
    if int >19
        
        if strfind(tline,'Axial')
            clips_coordAX = 1;
            clips_coordLat = 0;
            continue
        end
        
        if strfind(tline,'Lateral')
            clips_coordAX = 0;
            clips_coordLat = 1;
            continue
        end
        
        if clips_coordAX
            lineTxt = strsplit(tline);
            axial(count_clips).xy(count_points,:) = ...
                [str2num(lineTxt{1,size(lineTxt,2)-2}),str2num(lineTxt{1,size(lineTxt,2)-1})];
            count_points = count_points + 1;
            if count_points > 18
                count_clips = count_clips+1;
                count_points = 1;
                if count_clips > num_of_clips
                    count_clips = 1;
                    clips_coordAX = 0;
                    continue
                end
            else
                continue
            end
                
        elseif clips_coordLat
            lineTxt = strsplit(tline);
            lat(count_clips).xy(count_points,:) = ...
                [str2num(lineTxt{1,size(lineTxt,2)-1}),str2num(lineTxt{1,size(lineTxt,2)-2})];
            count_points = count_points + 1;
            if count_points > 18
                count_clips = count_clips+1;
                count_points = 1;
                if count_clips > num_of_clips
                    count_clips = 1;
                    clips_coordLat = 0;
                    continue
                end
            else
                continue
            end
        end
            
            
        
            
            
            
    end
        
end
fclose(fid);
   




end

