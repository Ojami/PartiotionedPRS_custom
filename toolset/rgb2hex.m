function hex = rgb2hex(rgb)
% a light version of rgb2hex https://se.mathworks.com/matlabcentral/fileexchange/46289-rgb2hex-and-hex2rgb

if max(rgb(:))<=1
    rgb = round(rgb*255); 
else
    rgb = round(rgb); 
end


hex(:,2:7) = reshape(sprintf('%02X',rgb.'),6,[]).'; 
hex(:,1) = '#';
hex = string(hex);

end