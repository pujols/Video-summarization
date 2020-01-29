function res = seqDPP_evaluate(videos_te, inds_te, num_order, approach_name, dataset)

cd ../Evaluation
kill_seqDPP(approach_name, inds_te, dataset);
make_folder_oracle(videos_te, approach_name, inds_te, dataset)
[O_S_detail, ~] = comparision_seqDPP({approach_name}, num_order, inds_te, dataset);
[~, True_RP, True_F1] = sample_seqdpp_comp(O_S_detail);
res = [True_F1, True_RP]; % F-score, Recall, Precision
cd ../codes
end