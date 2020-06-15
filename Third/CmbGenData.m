function [CmbData MidNum] = CmbGenData(GenData,M)

Num = size(GenData,1);
  
index = unique(GenData,'rows');
MidNum = (3^M+1)/2; 

for i=1:Num
    for ic = 1:length(index)

        if GenData(i,:) == index(ic,:)
            CmbData(i,:) = ic-MidNum;
            break
        end
    end
end        