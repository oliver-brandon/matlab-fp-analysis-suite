filename = '/Users/brandon/My Drive/prl/PRL_GRABDA/newCohortMats (Processed)/201F_Acq1_JZL8.mat';
load(filename)

jsonStr = jsonencode(data,'PrettyPrint',true);
jsonFileName = '201F_Acq1_JZL8.json';
fid = fopen(jsonFileName, 'w');
fprintf(fid, '%s', jsonStr);
fclose(fid);
