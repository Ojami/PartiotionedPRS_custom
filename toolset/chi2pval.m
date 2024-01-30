function p = chi2pval(x, nu)
% chi2 to p-value when p-value is close to 0 (survival function). 
% this is the same as 1 - chi2cdf(x, nu) but more accurate for extreme
% values
% this is equivalent to chi2.sf: https://docs.scipy.org/doc/scipy-0.7.x/reference/generated/scipy.stats.chi2.html
 
p =  gammainc(x/2, nu/2, 'upper');
end