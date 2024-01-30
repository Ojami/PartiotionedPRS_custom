function name = getRandomName(name, len, maxloop)
% generates a random string from an input string pattern name.
% Oveis Jamialahmadi, Sahlgrenska Academy, Oct 2021.

arguments
    name {mustBeTextScalar} % must be in format of NAME.EXTENSION, e.g. myfile.mat
    len {mustBeGreaterThan(len, 2)} = 15 % string length
    maxloop {mustBeGreaterThan(maxloop, 2)} = 10 % max number to check 
end

[~, name, fileExt] = fileparts(name);
name = name + ".";
charSet = ['a':'z',upper('a':'z'),'0':'9'];
randStr = string(charSet(randi(numel(charSet),1, len)));
cnt = 1;
while exist(name + randStr + fileExt, 'file')
    randStr = string(charSet(randi(numel(charSet),1, len)));
    cnt = cnt + 1;
    if cnt > maxloop
        error('cannot generate a unique random string :(')
    end
end

name = name + randStr + fileExt;

end