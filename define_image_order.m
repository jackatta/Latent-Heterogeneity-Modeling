% Define what biological conditions the images correspond to for each
% experiment. Must be done manually!
% 
% This information is used to label each image based on the experiment,
% cell line and treatment condition. rep_str function creates as many
% labels as technical replicates of a given condition in an experiment.

function [label_table] = define_image_order()

load_dir = 'images\TMR channel';
file_list = dir( fullfile(load_dir,'*.tif') ); % image already manually ordered

%% MANUAL PART
name_matrix{1} = rep_str({'DEX-','VEH','BASELINE','EIPA40','EIPA80'},... % conditions in experiment
                     [3,    5,      5,          5,      5]); % number technical replicates for that condition
name_matrix{2} = rep_str({'DEX-','EGF-','BASELINE','VEH','EIPA40','EIPA80'},...
                     [3,    6,      6,          6,      6,      5]) ;
name_matrix{3} = rep_str({'DEX-','EGF-','BASELINE','VEH','EIPA40','EIPA80'},...
                     [3,    5,      6,          6,      6,      6]) ; 
name_matrix{4} = rep_str({'DEX-','VEH','EGF-','BASELINE','EIPA40','EIPA80'},...
                     [3,    6,      5,          4,      5,      5]) ; 
for jj=5:9
    name_matrix{jj} = rep_str({'DEX-','EGF-','BASELINE','VEH','EIPA40','EIPA80'},...
                     [3,    6,      6,          6,      6,      6]) ;  
end

%% back to auto
for ii=1:length(name_matrix) % no. experiments
    for jj=1:length(name_matrix{ii}) % no. conditions
        reps = size(name_matrix{ii}{jj},1); % no. tech reps
        name_matrix{ii}{jj}(:,2) = repmat({file_list(ii).name},reps,1);
        name_matrix{ii}{jj}(:,3) = repmat(regexp(file_list(ii).name, '(Bx[PCGR3680]+)', 'match'),reps,1);
        name_matrix{ii}{jj}(:,4) = repmat(regexp(file_list(ii).name, '(5JH\d+)', 'match'),reps,1);
        name_matrix{ii}{jj}(:,5) = num2cell([1:reps]');
    end
    name_matrix{ii} = cat(1,name_matrix{ii}{:});
end
name_matrix = cat(1,name_matrix{:});
name_matrix = [name_matrix(:,2:4) name_matrix(:,1) name_matrix(:,5)] ; % file name, cell line, experiment, condition, tech rep no.

label_table = cell2table(name_matrix,'VariableNames',{'FileName','CellLine','Experiment','Treatment','TechRep'});
label_table.CellLine( strcmp(label_table.CellLine,'BxGR80C') ) = repmat({'BxGR080C'},sum(strcmp(label_table.CellLine,'BxGR80C')),1); % easier ordering with BxGR360C

end


%% helper function
function [cell_string] = rep_str(input_strings, n_vec)
if numel(input_strings)~=numel(n_vec)
    error('Cell array and N vector not same size');
end
cell_string = cell(1,length(n_vec));
for ii=1:length(n_vec)
%     cell_string{ii} = repmat({input_strings{ii}},1,n_vec(ii));
    cell_string{ii} = repmat({input_strings{ii}},n_vec(ii),1);
end
end
