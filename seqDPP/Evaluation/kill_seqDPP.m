function kill_seqDPP(approach_name, te_inds, dataset)

foldername = '../data/OVP_YouTube_cmp';

if (strcmp(dataset, 'OVP'))
    OVP_YouTube_index = 21 : 70;
elseif (strcmp(dataset, 'YouTube'))
    OVP_YouTube_index = [11 : 20, 71 : 110];
end

for n = 1 : length(te_inds)
    videoName = sprintf('v%d', OVP_YouTube_index(te_inds(n)));
    if (isdir(fullfile(fullfile(foldername, videoName), approach_name)))
        rmdir(fullfile(fullfile(foldername, videoName), approach_name), 's');
    end
end

end