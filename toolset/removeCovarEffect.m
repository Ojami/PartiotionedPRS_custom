function [X, y, ZSZX, ZSZy] = removeCovarEffect(X, y, Z)
% Taken from SuSiE remove.covariate.effects
% Oveis Jamialahmadi, University of Gothenburg. Dec 2022.

% @22NOV2023: now can handle categorical covariates (can be left as string
% in the input Z table) by converting them to dummy variables.

arguments
    X {mustBeNonempty,mustBeA(X, ["table", "double"])} % predictors
    y {mustBeNonempty, mustBeA(y, ["table", "double"])} % response
    Z {mustBeNonempty, mustBeA(Z, ["table", "double"])} % covariates to be removed
end

Xnames = [];
if istable(X)
    Xnames = X.Properties.VariableNames;
    X = X{:, :};
end

ynames = [];
if istable(y)
    ynames = y.Properties.VariableNames;
    y = y{:, 1};
end

% check if there are categorical covariates and convert them to dummy
% variables
if ~istable(Z)
    Z = array2table(Z); 
end

dt = varfun(@class, Z, "OutputFormat", "cell");
dt = string(dt);
strCovars = dt == "string";
if any(strCovars) % convert strings to categorical
    Z = convertvars(Z, strCovars, @categorical);
end

catCovars = ismember(dt, ["string", "categorical"]);
if any(catCovars)
    nc = find(catCovars);
    catcovar = Z(:, nc);
    Z(:, nc) = [];
    dcvar = cell(numel(nc), 1);
    for k = 1:numel(nc)
        dcvar{k} = dummyvar(catcovar{:, k});
        dcvar{k}(:, 1) = [];
    end
    Z = [Z{:, :}, horzcat(dcvar{:})];
end

if istable(Z), Z = Z{:, :}; end


if any(isnan(X), "all")
    error("removeCovarEffect: NaN values in X matrix are not allowed!")
elseif any(isnan(y), "all")
    error("removeCovarEffect: NaN values in y vector are not allowed!")
elseif any(isnan(Z), "all")
    error("removeCovarEffect: NaN values in Z matrix are not allowed!")
end

if (any(Z(:, 1) ~= 1))
    Z = [ones(size(Z, 1), 1), Z];
end

A = Z'*Z;
SZy = A\(y'*Z).';
SZX = A\(Z'*X);
ZSZy = Z*SZy;
ZSZX = Z*SZX;
y = y - ZSZy;
X = X - ZSZX;

if ~isempty(Xnames)
    X = array2table(X, VariableNames=Xnames);
end

if ~isempty(ynames)
    y = array2table(y, VariableNames=ynames);
end

end % END