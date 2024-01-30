function valin = irnt(valin)
% perform Inverse Normal Transformation (IRNT)
% Oveis Jamialahmadi, Sep. 2020. GU.

if nargin < 1
    error('no input?')
elseif ~isrow(valin) && ~isvector(valin)
    error('row/column vector only!')
end

% check nomrality 
if ~kstest(valin)
    warning('input data are already normal!')
end
valin = norminv((tiedrank(valin)-0.5)./sum(~isnan(valin)));
end % END