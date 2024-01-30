function [hdr, changedIdx] = createGWASheader(hdr, opts)
% creates a standard GWAS header from GWAS summary stats. This header can
% replace the hdr of an exsiting summary stat data to be used in downstream
% tools, e.g. PolyFun, FUMA or MAGMA.
% Oveis Jamialahmadi, University of Gothenburg, March 2023

arguments
    hdr {mustBeText, mustBeVector, mustBeNonmissing}
    opts.method {mustBeMember(opts.method, ["GWAS", "TwoSampleMR"])} = "GWAS" % table columns format 
end

hdr = string(hdr);
hdr_raw = hdr;

pcolidx = find(contains(hdr.lower, ...
    ["p_bolt_","p_value", "pval", "p_adj", "p-adj", "log10p"]) | hdr.lower == "p", 1);
posidx = find(ismember(hdr.lower, ["bp", "position", "pos", "base_pair_location"]) | contains(hdr.lower, "genpos"), 1);
betaidx = find(ismember(hdr.lower, ["beta", "b", "beta_y1", "beta.y1", "lnor"]), 1);
chridx = find(startsWith(hdr.lower, "chr"), 1);
afreqidx = find(ismember(hdr.lower, ["af_allele2", "a1freq", "a1_freq", "maf", "a2freq"]), 1);
seidx = find(startsWith(hdr.lower, ["se", "standard_error"]), 1);
snpidx = find(ismember(hdr.lower, ["snp", "varid", "rsid", "rsids", "id", "variant_id"]), 1);
a2idx = find(ismember(hdr.lower,  ["allele0", "a0", "allele2", "a2", "other_allele"]), 1);
a1idx = find(ismember(hdr.lower,  ["allele1", "a1", "effect_allele", "effectall"]), 1);
nidx = find(ismember(hdr.lower, "n"), 1);
ncaseidx = find(ismember(hdr.lower, "n.case"), 1);

if opts.method == "GWAS"
    hdr(pcolidx) = "P";
    hdr(posidx) = "BP";
    hdr(chridx) = "CHR";
    hdr(betaidx) = "BETA";
    hdr(afreqidx) = "A1FREQ";
    hdr(seidx) = "SE";
    hdr(a1idx) = "A1";
    hdr(a2idx) = "A2";
    hdr(nidx) = "N";
elseif opts.method == "TwoSampleMR"
    hdr(pcolidx) = "pval";
    hdr(posidx) = "position";
    hdr(chridx) = "chr";
    hdr(betaidx) = "beta";
    hdr(afreqidx) = "eaf";
    hdr(seidx) = "se";
    hdr(a2idx) = "other_allele";
    hdr(a1idx) = "effect_allele";
    hdr(nidx) = "samplesize";
    hdr(ncaseidx) = "ncase";
end

hdr(snpidx) = "SNP";
changedIdx = (hdr ~= hdr_raw);
end