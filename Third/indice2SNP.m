function snpInfo = indice2SNP(vec)
snp = importdata('/Users/yiw/Desktop/data/SNP.dat');
snps = regexp(snp{1}, '\s+', 'split');
len = length(snps);
snpInfo = cell(size(vec));
for i = 1:size(vec)
    snpInfo(i) = snps(vec(i));
end
        
    