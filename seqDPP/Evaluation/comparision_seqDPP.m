function [output_record, output_summary] = comparision_seqDPP(approach_Name, order_cmp, te_inds, dataset)
%% Setup
write_seqDPP_input(['input_wholeCmp' num2str(order_cmp) '.txt'], approach_Name, te_inds, dataset);
threshold = 0.5;
user_Index = [1 2 3 4 5];
if (strcmp(dataset, 'OVP'))
    OVP_YouTube_index = 21 : 70;
elseif (strcmp(dataset, 'YouTube'))
    OVP_YouTube_index = [11 : 20, 71 : 110];
end
OVP_YouTube_index = OVP_YouTube_index(te_inds);

%% VSUMM Eval
% write_VSUMM_input_wholeCmp('input_wholeCmp.txt');
disp(['java -jar CUS' num2str(order_cmp) '.jar -i input_wholeCmp' num2str(order_cmp) '.txt -o output_wholeCmp' num2str(order_cmp) '.txt -u ' num2str(length(user_Index))...
    ' -a ' num2str(length(approach_Name)) ' -t ' num2str(threshold)]);
system(['/usr/lib/jvm/jre-1.6.0/bin/java -jar CUS' num2str(order_cmp) '.jar -i input_wholeCmp' num2str(order_cmp) '.txt -o output_wholeCmp' num2str(order_cmp)...
    '.txt -u ' num2str(length(user_Index)) ' -a ' num2str(length(approach_Name)) ' -t ' num2str(threshold)]);
system(['chmod 777 input_wholeCmp' num2str(order_cmp) '.txt -R']);
system(['chmod 777 output_wholeCmp' num2str(order_cmp) '.txt -R']);
[output_record, output_summary] = read_seqDPP_output(['output_wholeCmp' num2str(order_cmp) '.txt'], OVP_YouTube_index, length(user_Index), length(approach_Name));
delete(['input_wholeCmp' num2str(order_cmp) '.txt']);
delete(['output_wholeCmp' num2str(order_cmp) '.txt']);

end