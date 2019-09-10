% Hypothesis tests on trends using Wald tests

function [p_table] = trend_hypothesis_testing(mdl,trends,opt)
%% Setup more complicated hypothesis tests using emm with cell line marginalized (easier to establish on this smaller table)
[emmTM] = emmeans(mdl,{'Treatment',opt.gmm.latent_var},'unbalanced','effects','asfit','Offset',opt.glme.emm_offset); 

% Does EGF stimulation have a trend-dependent effect on MPC?
% % (EGF+) - (EGF-) = 0 
L{1} = get_contrast(emmTM.table.Treatment,{'BASELINE','EGF-'},[1 -1],trends);
h0names{1} = 'egf';
hstring{1} = 'H0: There is no effect on MPC by EGF stimulation';

% Does veh ctrl (DMSO) have a trend-dependent effect on MPC?
% % (EGF+) - (VEH) = 0 
L{end+1} = get_contrast(emmTM.table.Treatment,{'BASELINE','VEH'},[1 -1],trends);
h0names{end+1} = 'veh';
hstring{end+1} = 'H0: There is no effect on MPC by DMSO treatment';

% Is there a trend-dependent effect on the linear dose-dependence of EIPA on MPC?
L{end+1} = get_contrast(emmTM.table.Treatment,{'VEH','EIPA40','EIPA80'},[1 0 -1],trends); % negative linear trend b/c linear dosing of inhibitor
h0names{end+1} = 'eipadose';
hstring{end+1} = 'H0: There is no effect of dose-dependent EIPA inhibition';

% Is there a trend-dependent effect of unstimulated MPC level? 
L{end+1} = get_contrast(emmTM.table.Treatment,{'EGF-'},[1],trends);
h0names{end+1} = 'unstim';
hstring{end+1} = 'H0: There is no effect on unstimulated levels of MPC';

%% hypothesis tests PER CELL LINE 
[emmCTM] = emmeans(mdl,{'CellLine','Treatment',opt.gmm.latent_var},'unbalanced','effects','asfit','Offset',opt.glme.emm_offset);
stnames = unique(emmCTM.table.CellLine,'stable');

Lec = cell(1,length(L));
for ii=1:length(stnames)
    for jj=1:length(L)
        x = find(L{jj}~=0); % use above contrasts
        Lc = zeros(1,size(L{jj},1)*height(emmCTM.table));
        Lc(x + (ii-1)*height(emmTM.table)) = L{jj}(x);
        Lc = reshape(Lc,size(L{jj},1),[]);
        if ~strcmp(h0names{jj},'unstim') % no contrasts yet defined for 'unstim', just extracting relevant coefficients so far
            H0.(h0names{jj}).(stnames{ii}) = contrasts_wald(mdl,emmCTM,Lc);
        end
        Lec{jj} = [Lec{jj}; Lc]; % keep each contrast
    end
end


%% hypothesis tests on RhoG-dependent effects
% these tests are contrasts of contrasts of contrasts
% assume even spacing of RhoG expression for trend tests
r_trends(:,strcmp(stnames,'BxPC3')) = -1; 
r_trends(:,strcmp(stnames,'BxGR080C')) = 0; 
r_trends(:,strcmp(stnames,'BxGR360C')) = 1; 
trendnames = {'linear'};

for rr=1:length(h0names)
    if ~strcmp(h0names{rr},'eipajoint') % contrast of joint tests isn't necessary
        for tt=1:size(r_trends,1)
            % hypothesis test
            Lt = sum(bsxfun(@times, Lec{rr}', r_trends(tt,:)),2)'; % trend contrast
            H0.(h0names{rr}).(trendnames{tt}) = contrasts_wald(mdl,emmCTM,Lt);
        end
    end
end

p_table = print_results(H0, hstring);

end

%% helper functions
function [L] = get_contrast(celstr,con_str,con_wgts,con_trends)

if length(con_str)~=length(con_wgts)
    error('con_str and con_wgts must be same dimensions.');
elseif sum(ismember(celstr,con_str{1}))~=length(con_trends)
    error('number of instances of a single con_str must be same as length of con_trends.');
elseif any(~cellfun(@(x) any(ismember(unique(celstr),x)), con_str))
    error('con_str probably spelled wrong. must be member of celstr.');
end

L = zeros(length(con_str),length(celstr));
for ii=1:length(con_str)
    L(ii,:) = double(strcmp(celstr,con_str{ii}))';
    L(ii,L(ii,:)~=0) = L(ii,L(ii,:)~=0) .* con_trends';
end
L = sum(bsxfun(@times,L',con_wgts),2)';

end

function [p_table] = print_results(H0,hstrings)
% get p-values of contrasts (but not trends)
hfields = fieldnames(H0);
cnames = fieldnames(H0.(hfields{1}));
p_array = cell(length(hfields),length(cnames));
for ii=1:length(hfields)
    for jj=1:length(cnames)
        try
            p_array{ii,jj} = H0.(hfields{ii}).(cnames{jj}).pVal;
        catch MERR
            if strcmp(MERR.identifier,'MATLAB:nonExistentField')
                p_array{ii,jj} = '--';
            else
                rethrow(MERR);
            end
        end
    end
end
p_table = cell2table(p_array, 'VariableNames', cnames', 'RowNames', hstrings);

end

