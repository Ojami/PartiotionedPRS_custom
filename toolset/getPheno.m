function pheno = getPheno(loadflag, mergeflag, pheno_dir)
% getPheno is a simple function
% able to search/get/load a set of target phenotypes from already saved
% phenotype mat files.
% OPTIONAL INPUT:
%   loadflag: whether load traits (true) or just return a string array
%   containing the name of traits (false, default).
%   mergeflag: whether to merge (ture) loaded traits into a single table
%   (default: false)
%   pheno_dir: phenotype folder. By fefault 'UKB_Pheno' folder in the same
%   directory as of getPheno function will be used.
% 
% OUTPUT:
%   pheno: a string array of trait names or a struct/table with loaded traits.
% 
% Oveis Jamialahmadi, Sahlgrenska Akademy, June 2021.
% 
% @15/09/2021: a bug fixed with 'mergeflag' to exclude 'exeid' prior to
%              merging.
% @20/09/2021: a bug was fixed with 'exeid' field.
% @15/10/2021: a bug was fixed with long tag names, when 'mergeflag' set to
%              true

if nargin < 2 || isempty(mergeflag)
    mergeflag = false;
end
if nargin < 1 || isempty(loadflag)
    loadflag = false;
end

if mergeflag
    loadflag = true;
end

if nargin < 3 || isempty(pheno_dir)
    phenofolder = 'UKB_PHENO'; % must be in path
    homeDir = which('getPheno.m');
    homeDir = fileparts(homeDir); 
    pheno_dir = fullfile(homeDir, phenofolder);
else
    mustBeFolder(pheno_dir)
end

pheno = {dir(fullfile(pheno_dir, '*.mat')).name}.';
pheno_names = regexprep(pheno, '.mat$', '');
if isempty(pheno_names)
    fprintf('getPheno: no phenotype was found in %s!\n', pheno_dir)
    return
end
sel_flag = true;
pheno_collector = [];
pheno_name_numbered = strcat(strsplit(num2str(1:numel(pheno_names))).','-', pheno_names);
while sel_flag % Loop over phenotype selection: possible to select several times
    fprintf('List of existing phenotypes:\n')
    fprintf('%s\n', pheno_name_numbered{:})
    pheno_sel = input('Insert indices or search by keyword (use ; to separate different phenotypes): ', 's');
    try
        pheno_collector = [pheno_collector, eval(pheno_sel)];
    catch
        if isempty(pheno_sel)
            pheno_collector = 1:numel(pheno_names);
            fprintf('All phenotypes were selected\n')
            break
        else % Character keywords
            pheno_sel = lower(strtrim(strsplit(pheno_sel, ';')));
            pheno_collector_inner = (0);
            for ii = 1:numel(pheno_sel) % Loop over input keywords
%                 found_pheno_idx = ismember(lower(pheno_names), pheno_sel(ii));
%                 if ~any(found_pheno_idx)
                    found_pheno_idx = contains(lower(pheno_names), pheno_sel(ii));
%                 end
                if sum(found_pheno_idx) > 1 % If matches > 1 pheno
                    found_phenos = pheno_names(found_pheno_idx);
                    fprintf('%s matches more than one phenotype:\n', pheno_sel{ii})
                    for jj = 1:numel(found_phenos)
                        fprintf('%d-%s\n', jj, found_phenos{jj})
                    end
                    found_phenos_sel = input('?[Enter to select all] ');
                    if isempty(found_phenos_sel)
                        found_phenos_sel = 1:numel(found_phenos);
                    end
                    if ~isempty(found_phenos_sel)
                        found_phenos_sel = found_phenos(found_phenos_sel);
                        selected_multiphenos = find(ismember(pheno_names, found_phenos_sel));
                        pheno_collector_inner(numel(pheno_collector_inner)+1:...
                            numel(pheno_collector_inner)+numel(selected_multiphenos)) =...
                            selected_multiphenos;
                    end

                elseif ~sum(found_pheno_idx) % If not matches was found
                else % If unique match is found
                    pheno_collector_inner(numel(pheno_collector_inner)+1) = find(found_pheno_idx);
                end
            end
            pheno_collector_inner(pheno_collector_inner < 1) = [];
            pheno_collector = [pheno_collector, pheno_collector_inner];
        end
    end
    pheno_continue = input('Do you want to select more phenotypes (y/[n])? ', 's'); 
    if ~strcmpi(pheno_continue, 'y')
        sel_flag = false;
    end

end
pheno_collector = unique(pheno_collector, 'stable');
pheno_names = pheno_names(pheno_collector);
pheno_names = strcat(strsplit(num2str(1:numel(pheno_names))).','-', pheno_names);
if isempty(pheno_names)
    fprintf('No phenotype was selected\n')
    return
end
fprintf('Selected phenotypes:\n')
fprintf('%s\n', pheno_names{:})

if loadflag
    pheno_names = fullfile(pheno_dir, pheno(pheno_collector));
    fprintf('\nLoading the selected phenotypes...\n')
    pheno = struct;
    if mergeflag; eid = []; end
    exeid = []; invalidIdx = [];
    for ii = 1:numel(pheno_names)
        load_pheno = load(pheno_names{ii});
        try
            pheno(ii).eid = load_pheno.UKB_STRUCT_ALL.eid;
        catch
            fprintf('ERROR::getPheno::%s is not a valid pheno struct!\n', pheno_names{ii})
            invalidIdx = [invalidIdx, ii];
            continue
        end
        if (load_pheno.UKB_STRUCT_ALL.numericFlag || strcmpi(load_pheno.UKB_STRUCT_ALL.tag, 'sex')) && mergeflag
            if isempty(eid)
                eid = pheno(ii).eid;
            else
                eid = intersect(eid, pheno(ii).eid);
            end
        end
        pheno(ii).rawUKB = load_pheno.UKB_STRUCT_ALL.rawUKB;
        pheno(ii).numericFlag = load_pheno.UKB_STRUCT_ALL.numericFlag;
        pheno(ii).termMeaning = load_pheno.UKB_STRUCT_ALL.termMeaning;
        pheno(ii).tag = load_pheno.UKB_STRUCT_ALL.tag;
        try
            pheno(ii).exeid = load_pheno.UKB_STRUCT_ALL.exeid;
            exeid = union(exeid, pheno(ii).exeid);
        catch
            pheno(ii).exeid = [];
        end
        
    end
    pheno(invalidIdx) = [];
    
    if mergeflag
        eid = setdiff(eid, exeid); % remove exeids 
        res = table(eid);
        for i = 1:numel(pheno)
            if pheno(i).numericFlag || strcmpi(pheno(i).tag, 'sex')
                [~, idx] = ismember(eid, pheno(i).eid);
                try
                    res.(pheno(i).tag) = double(pheno(i).rawUKB(idx));
                catch
                    trname = char(matlab.lang.makeValidName(pheno(i).tag));
                    res.(trname) = double(pheno(i).rawUKB(idx));
                end
            else
                idx = ismember(eid, pheno(i).eid);
                try
                    res.(pheno(i).tag) = zeros(numel(eid), 1);
                    res.(pheno(i).tag)(idx) = 1;
                catch % long tag names may exceed the max allowed length for variable names
                    trname = char(matlab.lang.makeValidName(pheno(i).tag));
                    res.(trname) = zeros(numel(eid), 1);
                    res.(trname)(idx) = 1;
                end
            end
        end
        pheno = res;
    end
    fprintf('\b\b Done!\n')
else
    pheno = pheno(pheno_collector);
    pheno = string(regexprep(pheno, '.mat$', ''));
end


end % END