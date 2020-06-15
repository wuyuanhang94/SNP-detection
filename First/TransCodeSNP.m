
close all
clear
clc

A = importdata('E:\gmcm\genotype.dat');
for ic = 1:1001
    S(ic,:) = regexp(A{ic}, '\s+', 'split');
end

for ir = 1:9445
    index = unique(S(2:end,ir));
    
    for ic = 1:1000
        switch S{ic+1,ir}
            case index{1}
                GenData{ic,ir} = '00';%-1;
            case index{2}
                GenData{ic,ir} = '01';
            case index{3}
                GenData{ic,ir} = '11';
        end
    end 
end
% save('GenData.mat','GenData')