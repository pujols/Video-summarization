function [output_record, output_summary] = read_seqDPP_output(file_name, video_Index, Num_user, Num_approach)

%% Input:
% file_name: The file name to read. ex: 'output.txt'
% video_Index: The videos that you want to extract the information from.
%              ex: (21:70) for all the videos
% Num_user: The number of users considered. ex: 5 (in general)
% Num_approach: The number of approaches compared


%% Data setting
output_record = cell(Num_approach, length(video_Index));
% Each cell element is a # users x 2 (CUSa, CUSe) matrix
output_summary = zeros(Num_approach, 2);
% Each row contains the (CUSa, CUSe) of that mehtod

%% file setting
fid = fopen(file_name);
now_line = fgetl(fid);

% while((length(now_line)~=1) || (sum(now_line==-1)~=1)) %% Not at the end of the document
while(isempty(strfind(now_line, '##########'))) %% Not at the end of the document
    if (~isempty(strfind(now_line, '## VideoName'))) %% New video starting

%% Locating the video: Method 1
%         video_important_check = zeros(1,length(video_Index));
%         for n = 1:length(video_Index) %% Checking video index
%             if (length(strfind(now_line, num2str(video_Index(n)))) == 1)
% %                 now_video = n;
% %                 break;
%                 video_important_check(n) = 1; 
%             end
%         end
%         video_Index_important_check = video_Index(video_important_check ~= 0);
%         now_video_important_check = max(video_Index_important_check);
%         now_video = find(video_Index == now_video_important_check);
%% Locating the video: Method 2
        video_name_important_loc1 = strfind(now_line, '/');
        video_name_important_loc2 = strfind(now_line, '\');
        if (isempty(video_name_important_loc1))
            video_name_important = str2double(now_line(max(video_name_important_loc2)+2:end));
        elseif (isempty(video_name_important_loc2))
            video_name_important = str2double(now_line(max(video_name_important_loc1)+2:end));
        else
            if (max(video_name_important_loc1) > max(video_name_important_loc2))
                video_name_important = str2double(now_line(max(video_name_important_loc1)+2:end));
            else
                video_name_important = str2double(now_line(max(video_name_important_loc2)+2:end));
            end
        end
        now_video = find(video_Index == video_name_important);
%% END        
        approach_record = zeros(Num_user,2);
        
    elseif (~isempty(strfind(now_line, 'UserSummary'))) %% Starting user summary
        now_user = str2double(now_line(length('UserSummary')+1));
%         loc_Approach = strfind(now_line, 'Approach');
%         now_Approach = str2double(now_line(length('Approach')+loc_Approach));
        loc_num = strfind(now_line, '.');
        if (now_line(loc_num(1)+2) ~= ' ')
            CUSa = str2double(now_line(loc_num(1)-1:loc_num(1)+2));
        else
            CUSa = str2double(now_line(loc_num(1)-1:loc_num(1)+1));
        end
        CUSe = str2double(now_line(loc_num(2)-1 : min(loc_num(2)+2,length(now_line))));
        approach_record(now_user,:) = [CUSa, CUSe];
        
    elseif (~isempty(strfind(now_line, 'Approach'))) %% End of an approach
        if (now_line(length('Approach')+2) ~= ' ')
            now_Approach = str2double(now_line(length('Approach')+1 : length('Approach')+2));
        else
            now_Approach = str2double(now_line(length('Approach')+1));
        end
        output_record{now_Approach, now_video} = approach_record;
    end
  
    now_line = fgetl(fid);
end

now_line = fgetl(fid);
for n = 1:Num_approach
    loc_num = strfind(now_line, '.');
    if (now_line(loc_num(1)+2) ~= ' ')
        CUSa = str2double(now_line(loc_num(1)-1:loc_num(1)+2));
    else
        CUSa = str2double(now_line(loc_num(1)-1:loc_num(1)+1));
    end
    CUSe = str2double(now_line(loc_num(2)-1 : min(loc_num(2)+2,length(now_line))));
    output_summary(n,:) = [CUSa, CUSe];
    now_line = fgetl(fid);
end

fclose(fid);

end