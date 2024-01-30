function [fid, msg] = runbash(cmd, bashname, opts)
% writes a input string vector, cmd, to a bash script, bashname, and runs
% it via wsl.

% Oveis Jamialahmadi, University of Gothenburg, March 2022.
% 
% @07/22/2022: 'conda' option was added for running the scripts under a
%               certain conda environment (conda info --envs to see all
%               envs).

arguments
    cmd {mustBeText, mustBeVector}
    bashname {mustBeTextScalar}
    opts.delete (1,1) logical = true % delete bash file after run
    opts.dir {mustBeFolder} = pwd % directory to save/run bash script to/from
    opts.waitprocess {mustBeTextScalar} = "" % process to wait for (e.g. R)
    opts.wait (1,1) logical = true % add wait to bash file 
    opts.parallel (1,1) logical = true % run all in the background
    opts.conda {mustBeTextScalar} % name of the conda environment
    opts.verbose (1,1) logical = false
    opts.runfrom {mustBeTextScalar} % run from this directory (WSL2 is usually faster is run from Linux directory); e.g. $HOME
end

opts.dir = string(opts.dir);
opts.dirwsl = makeWSLpath(opts.dir);
if ~endsWith(opts.dirwsl, "/"); opts.dirwsl = opts.dirwsl + "/"; end

cmd = string(cmd);
if isrow(cmd); cmd = cmd'; end
cmd(cmd == "" | ismissing(cmd)) = [];

if isfield(opts, 'conda') && opts.conda ~= ""
    cmd = ["conda activate " + opts.conda; cmd];
    opts.conda = true;
else
    opts.conda = false;
end

if opts.parallel
    cmd = regexprep(cmd, "&$", "");
    cmd = cmd + " &";
end

if opts.wait; cmd = [cmd; "wait"]; end % wait for all commands to finish

bashname = string(bashname);
if ~endsWith(bashname, '.sh'); bashname = bashname + ".sh"; end
writematrix(["#!/bin/bash"; cmd], fullfile(opts.dir, bashname), 'FileType', 'text', 'QuoteStrings', false);

if ispc
    cmdpref = "wsl ";
    [~, ~] = system(cmdpref + 'dos2unix "' + opts.dirwsl + bashname + '"');
else
    cmdpref = "";
end

if isfield(opts, "runfrom")
    if ~endsWith(opts.runfrom, "/"), opts.runfrom = opts.runfrom + "/"; end
    system(cmdpref + 'mv "' + opts.dirwsl + bashname + '" "' + opts.runfrom + bashname + '"');
    opts.dirwsl = opts.runfrom;

end

if ispc
    % check network drives mount, if they're not mounted re-apply mounting
    [~, msg] = system('net use');
    msg = splitlines(string(msg));
    msg(~startsWith(msg, "OK")) = [];

    % keep only mounted drives
    idx = contains(msg, "OK" + whitespacePattern + lettersPattern + ":");
    msg(~idx) = [];
    
    net.drive = extractBetween(msg, "OK" + whitespacePattern, ":" + whitespacePattern);
    net.path = extractBetween(msg, ":" + whitespacePattern, whitespacePattern);
    for i = 1:numel(net.drive)
        [~, msg] = system("wsl ls /mnt/" + lower(net.drive(i)));
        if isempty(msg) || string(msg).contains("No such file or directory")
            [~, ~] = system("wsl sudo mkdir -p /mnt/" + lower(net.drive(i)));
            [~, ~] = system('wsl sudo mount -t drvfs "\' + net.path(i) + '" /mnt/' + lower(net.drive(i)));
        end
    end
end

if opts.conda
    [fid, msg] = system(cmdpref + '/bin/bash -i "' + opts.dirwsl + bashname + '"');
    if fid
        [~, msg] = system(cmdpref + 'bash -i "' + opts.dirwsl + bashname + '"');
    end
else
    [fid, msg] = system(cmdpref + '/bin/bash -ilc "' + opts.dirwsl + bashname + '"');
    if fid
        [fid, msg] = system(cmdpref + 'bash -ilc "' + opts.dirwsl + bashname + '"');
        if fid
            [~, msg] = system(cmdpref + '"' + opts.dirwsl + bashname + '"');
        end
    end
end

if opts.verbose, disp(msg); end
% system("wsl systemd-run --scope -p MemoryLimit=1000M " + opts.dirwsl + bashname);

if opts.waitprocess ~= "" % wait for the commands to be done
    pause(5)
    [~, isRunning] = system(cmdpref + "pgrep " + opts.wait + " && wsl echo Running");
    while ~isempty(isRunning)
        [~, isRunning] = system(cmdpref + "pgrep -x " + opts.wait + " && wsl echo Running");
    end
end

if opts.delete
    if isfile(fullfile(opts.dir, bashname)), delete(fullfile(opts.dir, bashname)); end
    if isfield(opts, "runfrom")
        system(cmdpref + 'rm "' + opts.dirwsl + bashname + '"');
    end
end

end 