clear all
%Cbrain Variables%
a = 9.65;
alpha = 0.097;
beta = 0.642;
EXCELPATH = '/Users/brandon/Desktop/Book1.xlsx';
SHEET = 1;
RANGE = "D3:D250";
exceldata = readcell(EXCELPATH, 'Sheet', SHEET, 'Range', RANGE);
infusion_ts = cell2mat(exceldata(3:end-2, 1));
end_ts = height(infusion_ts);

K = 30;

for i = 1:end_ts
    if i == end_ts
        break
    end
    t = infusion_ts(i,1);
    u = infusion_ts(i+1,1);
    ts_new = double(t:K:u).';
    Cbrain(1:numel(ts_new),i) = ts_new;
end
Cbrain = Cbrain/60;
Cbrain = infusion_ts/60;
length_Cbrain = length(Cbrain);


for j = 1:length_Cbrain
    if j == length_Cbrain
        break
    end
    for m = 2:height(Cbrain)
        if m == height(Cbrain)
            j = j + 1;
        end
        m1 = Cbrain(m,j)-Cbrain(1,j);
        ITI(m,j) = m1; 
        
    end
end
for z = 3:length(ITI)
    z1 = Cbrain(1,z) - Cbrain(1,z-1);
    ITI(1,z) = z1;
end
Csum = ITI(ITI>=0);    
end_ts2 = height(Csum);
for y = 1:end_ts2
    diff1 = Csum(y,1);
    XX = a*(exp((-alpha)*diff1)-exp((-beta)*diff1));
    Csum(y,2) = XX;
end
cocaine = 2.192808952246444;
[r, index] = min(abs(Csum(1:end,2) - cocaine));