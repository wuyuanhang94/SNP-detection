importance = importdata('/Users/yiw/Desktop/data/importance_RF.mat');
genData = importdata('/Users/yiw/Desktop/data/GenData.mat');
[y,indices]=sort(importance,'descend');
subImpor = [indices(1:200,1),y(1:200,1)];
weights = [subImpor(:,1),subImpor(:,2)/sum(subImpor(:,2))];
weight = zeros(300,1);

for count = 1:300
    w = 0;
    j = 0;
    path = ['/Users/yiw/Desktop/trans_gene_info/gene_',num2str(count),'.mat'];
    gene = importdata(path);
    len = length(gene);
    for i = 1:len
        index = find(weights(:,1) == gene(i));
        if isempty(index)
            w = w;
        else
            j = j+1;
            %w = w + weights(index,2);
            buffer(j,1) = gene(i);
            buffer(j,2) = weights(index,2);
        end
    end
    if j>0
            [~,ibuffer] = max(buffer(:,2));
            imax = buffer(ibuffer,1);
            maxImpor = buffer(ibuffer,2);
            buffer(ibuffer,:) = [];
            other = buffer;
    end
    if j>1
        for k = 1:length(other(:,1))
            maxsnp = genData(:,imax);
            othersnp = genData(:,other(k,1));
            covValue = cov(maxsnp,othersnp);
            co = abs(covValue(1,1));
            other(k,2) = other(k,2)*(1-co);
        end
        weight(count) = (sum(other(:,2))+maxImpor);
    end
    if j == 1
        weight(count) = maxImpor;
    end
    if j == 0
        weight(count) = 0;
    end 
end