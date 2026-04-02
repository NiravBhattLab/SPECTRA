clear
% the consistent universal AGORA model
load('..\gapFilledAGORA_DB\ConsUmodel.mat')

p = dir('results/inputs');
p = {p(3:end).name}';

model_names = cellfun(@(x) strrep(x, '_maxsep.csv', ''), p, 'UniformOutput', false);

% empty table to store the results
exact_match_table = table('Size', [numel(model_names), 3], 'VariableTypes', {'string', 'double', 'double'}, 'VariableNames', {'Model', 'No_of_new_rxns', 'No_of_exact_match'});
for m = 1:numel(model_names)
    model_name = model_names{m};   
    [n_new_rxns,no_of_exact_match] = FindAccOfNewReactions(ConsUmodel,model_name);
    exact_match_table.Model{m} = model_name;
    exact_match_table.No_of_new_rxns(m) = n_new_rxns;
    exact_match_table.No_of_exact_match(m) = no_of_exact_match;
end

writetable(exact_match_table,'New_AGORA_clean_Results.csv')
function [n_new_rxns,no_of_exact_match] = FindAccOfNewReactions(ConsUmodel,model_name)

    % loading the new AGORA2 model
    new_model = load(['..\gapFilledAGORA_DB\AGORA_new_db\',model_name,'.mat']);
    new_model = new_model.model;


    % loading the exisiting AGORA2 model
    model = load(['..\gapFilledAGORA_DB\AGORA2\',model_name,'.mat']);
    model = model.model;


    new_rxns = setdiff(new_model.rxns,model.rxns);
    n_new_rxns = numel(new_rxns);

    ec_nos = ConsUmodel.rxnECNumbers(ismember(ConsUmodel.rxns,new_rxns));
    % removing empty entries
    ec_nos = ec_nos(~cellfun('isempty',ec_nos));


    dl_pred_ec = readtable(['results/inputs/',model_name,'_maxsep.csv'], 'Delimiter', ',','ReadVariableNames', false);

    % getting the column names
    col_names = dl_pred_ec.Properties.VariableNames;
    pred_ec ={};
    for c =2:length(col_names)
        col = col_names{c};
        col = dl_pred_ec.(col);
        for k=1:numel(col)
            ec_temp = split(col{k}, '/');
            ec_temp = ec_temp{1};
            if contains(ec_temp, 'EC')
                pred_ec = [pred_ec; strrep(ec_temp, 'EC:', '')];
            end
        end
    end

    no_of_exact_match = 0;

    for i=1:length(ec_nos)
        ec_temp = split(ec_nos{i}, ',');
        for k=1:length(ec_temp)
            if ismember(ec_temp{k}, pred_ec)
                no_of_exact_match = no_of_exact_match + 1;
                break
            end
        end
    end


    % pre_ec_sub = split(pred_ec,'.');
    % pre_ec_sub_temp = {};
    % 
    % for i=1:size(pre_ec_sub,1)
    %     pre_ec_sub_temp = [pre_ec_sub_temp; strjoin(pre_ec_sub(i,1:3),'.')];
    % end
    % 
    % 
    % no_of_sub_match = 0;
    % for i=1:length(ec_nos)
    %     ec_temp = split(ec_nos{i}, ',');
    %     for k=1:length(ec_temp)
    %         ec_temp_sub = split(ec_temp{k}, '.');
    %         ec_temp_sub = strjoin(ec_temp_sub(1:3),'.');
    %         if ismember(ec_temp_sub, pre_ec_sub_temp)
    %             no_of_sub_match = no_of_sub_match + 1;
    %             break
    %         end
    %     end
    % end
end


