function pth2 = fileparts2(pth, upn)
% similar to fileparts but returns the upper directory based on upn: upper
% directories
% Oveis Jamialahmadi, University of Gothenburg, August 2023

arguments
    pth {mustBeFolder}
    upn (1,1) double % number of upper directories 
end

pth2 = pth;
if upn <= 0 
    return
end

for k = 1:upn
    pth2 = fileparts(pth2);
end

end % END