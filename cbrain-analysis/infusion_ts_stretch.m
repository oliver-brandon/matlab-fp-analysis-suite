EXCELPATH = '/Users/brandon/My Drive (boliv018@ucr.edu)/DA WIN/VEH SA - CocBrain.xlsx';
SHEET = 7
RANGE = "AN1:AN14"
exceldata = readcell(EXCELPATH, 'Sheet', SHEET, 'Range', RANGE);
infusion_ts = cell2mat(exceldata(3:end-2, 1));
end_ts = height(infusion_ts);
x = 30; %time bins (s) to stretch the infusion timstamps%
j = 1;
ts_final = [];
for i = 1:end_ts
    t = infusion_ts(i,1);
    u = infusion_ts(i+1:end,1);
    ts_new = double(t:x:u).';
    len = height(ts_new);
    ts_stretch(j:j+len-1) = ts_new;
    ts_final = ts_stretch.';
    j = j+len;

end

for i = 1:end_ts
    t = infusion_ts(i,1);
    u = t+30;
    if u < infusion_ts(i+1,1)
        ts_final(i,1) = t;
        ts_final(i,2) = u;
        if i == end_ts-1
            break
        end
    else
        ts_final(i,1) = t;
    end
    
   
end

ts_final = cat(1,ts_final(:,1),ts_final(:,2));
ts_final = sort(ts_final);

        
    
    