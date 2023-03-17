clear all
%Cbrain Variables%
a = 9.65;
alpha = 0.097;
beta = 0.642;
EXCELPATH = '/Users/brandon/Desktop/Book1.xlsx';
SHEET = 1;
RANGE = "A1:A14";
exceldata = readcell(EXCELPATH, 'Sheet', SHEET, 'Range', RANGE);
infusion_ts = cell2mat(exceldata(3:end-2, 1));
end_ts = height(infusion_ts);

Cbrain = infusion_ts/60;
length_Cbrain = length(Cbrain);


for j = 2:length_Cbrain
    if j == length_Cbrain
        break
    end
    m = Cbrain(j,:);
    m1 = Cbrain(j+1,:);
    Csum(j) = (m1-m); 
    
end
Csum = Csum.';
for y = 1:height(Csum)
    diff1 = Csum(y,1);
    XX = a*(exp((-alpha)*diff1)-exp((-beta)*diff1));
    if y > 1
        XX = XX + Csum(y-1,2);
    end
    Csum(y,2) = XX;
end
