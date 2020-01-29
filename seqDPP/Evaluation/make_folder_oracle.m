function make_folder_oracle(videos, approach_name, te_inds, dataset)

foldername = '../data/OVP_YouTube_cmp';
framePath = '../data/Frames_sampled';


if (strcmp(dataset, 'OVP'))
    load ../video_summarization/Oracle_groundset/Oracle_OVP.mat Oracle_record
    OVP_YouTube_index = 21 : 70;
elseif (strcmp(dataset, 'YouTube'))
    load ../video_summarization/Oracle_groundset/Oracle_Youtube.mat Oracle_record
    OVP_YouTube_index = [11 : 20, 71 : 110];
end

for n = 1:length(te_inds)
    videoName = sprintf('v%d', OVP_YouTube_index(te_inds(n)));
    load(fullfile(framePath, [videoName, '.mat']));
    inds_frame = Oracle_record{te_inds(n), 3};
    inds_frame = intersect(15*(1:vidFrame.nrFramesTotal), sort(inds_frame));
    
    mkdir(fullfile(fullfile(foldername, videoName), approach_name));
    directPath = fullfile(fullfile(foldername, videoName), approach_name);
    
    frame_index = inds_frame(videos(n).Ypred);
    for m = 1:length(frame_index)
        imwrite(vidFrame.frames(frame_index(m)/15).cdata, fullfile(directPath, ['Frame' num2str(frame_index(m)) '.jpeg']), 'jpeg');
    end
    
    system(['chmod 777 ' directPath ' -R']);
end
end