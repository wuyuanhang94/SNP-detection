Genodata = importdata('/Users/yiw/Desktop/data/GenData.mat');

genChiMean = zeros(300,1);
genChiMax = zeros(300,1);
M = 3;
for count = 1:300
    path = ['/Users/yiw/Desktop/trans_gene_info/gene_',num2str(count),'.mat'];
    gene = importdata(path);
    len = length(gene);
    perm = combntns(gene,M);
    fitbuffer = zeros(size(perm,1),1);
    for i = 1:size(perm,1)
        cmb = Genodata(:,perm(i,:));
        [CmbData,MidNum] = CmbGenData(cmb,M);
        fitvalue = fitFunChi2(CmbData,MidNum);
        fitbuffer(i) = fitvalue;
    end 
    genChiMean(count) = mean(fitbuffer);
    genChiMax(count) = max(fitbuffer);
    
    clear fitbuffer
    clear gene
    clear cmb
    clear CmbData
end
    
            


