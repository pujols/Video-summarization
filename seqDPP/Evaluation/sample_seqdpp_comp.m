function [True_CU, True_RP, True_F1] = sample_seqdpp_comp(O_S_detail)

NUM_2 = size(O_S_detail, 1); % # approaches
True_CU = zeros(NUM_2, 2);
True_RP = zeros(NUM_2, 2);
True_F1 = zeros(NUM_2, 1);
   
for num2 = 1 : NUM_2
    now_OS2 = O_S_detail(num2, :);
    for now_at = 1:length(now_OS2) % for each video
        now_user_CU = now_OS2{now_at};
        per_video_CU = sum(now_user_CU, 1)/size(now_user_CU, 1);
        per_video_true_p = sum(now_user_CU(:,1)./(sum(now_user_CU, 2)+10^-10), 1)/size(now_user_CU, 1);
        per_video_true_F1 = sum((2*now_user_CU(:,1).*(now_user_CU(:,1)./(sum(now_user_CU, 2)+10^-10))) ./...
            (now_user_CU(:,1) + (now_user_CU(:,1)./(sum(now_user_CU, 2)+10^-10)) + 10^-10), 1)/size(now_user_CU, 1);
        True_CU(num2, :) = True_CU(num2, :) + per_video_CU/length(now_OS2);
        True_RP(num2, :) = True_RP(num2, :) + [per_video_CU(1), per_video_true_p]/length(now_OS2);
        True_F1(num2) = True_F1(num2) + per_video_true_F1/length(now_OS2);
    end
end


end