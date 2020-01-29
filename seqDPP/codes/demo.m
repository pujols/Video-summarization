%% Run seqDPP-linear on OVP (40 videos for training and 10 for testing) and YouTube (31 videos for training and 8 for testing)
rng_seed = 100;
display('#### seqDPP-linear on OVP ####');
seqDPP_linear_all('OVP', rng_seed);
display('Press any key to continue for seqDPP-linear on YouTube');
pause();

display('#### seqDPP-linear on YouTube ####');
seqDPP_linear_all('YouTube', rng_seed);
display('Press any key to continue for seqDPP-NN on OVP');
pause();

%% Run seqDPP-NN on OVP (40 videos for training and 10 for testing) and YouTube (31 videos for training and 8 for testing)
display('#### seqDPP-NN on OVP ####');
seqDPP_NN_all('OVP', rng_seed);
display('Press any key to continue for seqDPP-NN on YouTube');
pause();

display('#### seqDPP-NN on YouTube ####');
seqDPP_NN_all('YouTube', rng_seed);
pause();