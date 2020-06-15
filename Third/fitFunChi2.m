function fitValue = fitFunChi2(ChrValue,MidNum)

Num = size(ChrValue,1);
IndNum = MidNum*2-1;

M = zeros(2,IndNum);
X = zeros(1,IndNum);
Y = zeros(2,1);

data0 = ChrValue(1:500,:);
data1 = ChrValue(501:Num,:);
    
for i = 1:IndNum
    M(1,i) = length(find(data0 == i-MidNum));  
    M(2,i) = length(find(data1 == i-MidNum));
end
    
for i=1:2
        for k=1:IndNum
            if M(i,k) == 0 
                M(i,k)=M(i,k)+1;
            end
        end
end
    
X= sum(M);
Y = sum(M,2);
tot = sum(X);
    
E = Y*X/tot; 
fitValue = sum(sum(((M-E).^2)./E));


    

