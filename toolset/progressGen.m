function progressBar = progressGen(N, counter)
% generates a text progress bar for N calculations, so that counter can
% be 1:N. 
% Usage example:
% N = 100;
% textBar = progressGen(N); % initialize bar generator
% for i = 1:N
%   % do something
%   progressGen(textBar, i)
% end
% 
% Oveis Jamialahmadi, Sahlgrenska Academy, June 2021.

if nargin < 1
    error('progressGen: wrong number of inputs!')
end

inFlag = false;
if nargin < 2
    if numel(N) > 1 % wrong initializing, maybe counter is missing and first arg is progressBar?
        error('progressBar: first argument must be a scalar!')
    elseif isempty(N) || ~isnumeric(N)
        error('progressBar: first argument must be a nonempty double!')
    else
        inFlag = true; % initialize the progress bar
    end
end

if inFlag
    progressBar = floor(linspace(1, N+1, 11));
    fprintf(' [           ]') % generate an empty progress bar 
else % update progress bar
    progressTxt = [repmat('=', 1, sum(counter >= N)),'>',...
         repmat(' ', 1, 10-sum(counter >= N))];
     fprintf(repmat('\b', 1, 12))
     fprintf('%s]', progressTxt)
end


end