SNP = 'rs10779763';
snp = importdata('/Users/yiw/Desktop/data/SNP.dat');
chi2 = importdata('/Users/yiw/Desktop/mathmatical modeling/matlab/singleSNP/Chi2SingleSNP.mat');
snps = regexp(snp{1}, '\s+', 'split');
len = length(snps);
for i = 1:len
    if(strcmp(snps{i},SNP))
        chi = chi2(i)
        break
    end
end





