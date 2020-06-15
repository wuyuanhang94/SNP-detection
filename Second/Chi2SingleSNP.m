GenData = importdata('/Users/yiw/Desktop/mathmatical modeling/GenData.mat');
[m n]=size(GenData);
chi = zeros(n,1);
for i=1:1:n
    data0=GenData(1:500,i);
    data1=GenData(501:m,i);
    x00 = sum(data0 == -1);
    x01 = sum(data0 == 0);
    x02 = sum(data0 == 1);
    
    x10 = sum(data1 == -1);
    x11 = sum(data1 == 0);
    x12 = sum(data1 == 1);
    
    x0 = x00 + x10;
    x1 = x01 + x11;
    x2 = x02 + x12;
    
    y0 = x00 + x01 + x02;
    y1 = x10 + x11 + x12;
    
    tot = y0 + y1;
    
    z00 = y0*x0/tot;
    z01 = y0*x1/tot;
    z02 = y0*x2/tot;
    
    z10 = y1*x0/tot;
    z11 = y1*x1/tot;
    z12 = y1*x2/tot;
    
    X = [x00,x01,x02;x10,x11,x12];
    Z = [z00,z01,z02;z10,z11,z12];
    L = ((X-Z).^2)./Z;
    c = sum(L(:));
    chi(i) = c;
end
    