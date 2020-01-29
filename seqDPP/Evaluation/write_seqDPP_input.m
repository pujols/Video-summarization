function [] = write_seqDPP_input(filename, approach_Name, te_inds, dataset)

if (strcmp(dataset, 'OVP'))
    OVP_YouTube_index = 21 : 70;
elseif (strcmp(dataset, 'YouTube'))
    OVP_YouTube_index = [11 : 20, 71 : 110];
end
user_Index = 1:5;
foldername = '../data/OVP_YouTube_cmp';

fid = fopen(filename, 'w');

for k = 1:length(te_inds)
    videoName = sprintf('v%d',OVP_YouTube_index(te_inds(k)));
    fprintf(fid,['%s' videoName '\n'],[foldername '/']);

    for n = 1:length(user_Index)
        userName = sprintf('user%d',user_Index(n));
        fprintf(fid,['%s' videoName '%c' userName '\n'], [foldername '/'], '/');
    end
    
    for n = 1:length(approach_Name)
        fprintf(fid,['%s' videoName '%c' approach_Name{n} '\n'], [foldername '/'], '/');
    end
    fprintf(fid,'\n');
end
fclose(fid);

system(['chmod 777 ' filename ' -R']);
end




















