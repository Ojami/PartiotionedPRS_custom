function rcode = addRfont(rcode, opts)
% uses 'extrafont' R package to change font name on R graphs. By defaul
% looks for png() and pdf() functions.

arguments
    rcode {mustBeText}
    opts.fontFamily {mustBeTextScalar} = "Garamond" % adds font family only if finds a line in the inpuct script exporting figure as png/pdf
end

idx = find(rcode.startsWith(["png", "pdf", "cairo_pdf", "grDevices::cairo_pdf"] + "("));
for k = 1:numel(idx)
    rcode(idx(k)) = regexprep(rcode(idx(k)), ")$", ", family='" + opts.fontFamily + "')");
end

ercode(1, 1) = 'if(!"' + opts.fontFamily + ...
    '" %in% extrafont::fonts()) ' +  ...
    'font_import(pattern="' + opts.fontFamily + '", prompt = F)';
ercode(2, 1) = 'extrafont::loadfonts(device = "all")';

if isrow(rcode)
    rcode = [ercode', rcode];
else
    rcode = [ercode; rcode];
end


end