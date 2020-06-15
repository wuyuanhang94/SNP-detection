%clear
%clc
%genotype = importdata('/Users/yiw/Desktop/mathmatical modeling/genotype.dat');
snps = regexp(genotype{1}, '\s+', 'split');
path = '/Users/yiw/Desktop/gene_info/';
for count = 1:300
    genex = ['gene_',num2str(count),'.dat'];
    genePath = [path,genex];
    gene = importdata(genePath);
    len = length(gene);
    trans = zeros(len,1);
    for i = 1:len
        for j = 1:length(snps)
            genei = gene{i};
            snpsj = snps{j};
            if strcmp(genei,snpsj)
               trans(i) = j;
               break;
            end
        end
    end
    save(['/Users/yiw/Desktop/trans_gene_info/gene_',num2str(count),'.mat'],'trans');
    clear gene
end