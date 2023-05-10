function [stat1,stat2,stat3,stat4] = NERD_STATS(toc,numFiles)

% filepath = '../mat-stats';
% if ~exist(filepath,'nerd-stats.txt')
%     fid = fopen(strcat(filepath,'/','nerd-stats.txt'),'w');
%     fprintf(fid,'cumulative files processed: 0\n');
%     fprintf(fid,'cumulative hours runtime: 0\n');
%     fclose(fid);
%     filesCumulative = 0;
%     runtimeCumulative = 0;
% else
    % Read the current value of the cumulative number of files processed
    fid = fopen('../mat-stats/nerd-stats.txt', 'r');
    if fid ~= -1
        str = fgets(fid);
        tokens = strsplit(str, ':');
        filesCumulative = str2double(tokens{2});
        str = fgets(fid);
        fclose(fid);
        tokens = strsplit(str, ':');
        runtimeCumulative = str2double(tokens{2});
    else
        filesCumulative = 0;
        runtimeCumulative = 0;
    end
% end
runHrs = (toc/60)/60;

filesCumulative = filesCumulative + numFiles;
runtimeCumulative = runtimeCumulative + runHrs;

% Write the updated value to the text file
fid = fopen('../mat-stats/nerd-stats.txt', 'w');
fprintf(fid, 'cumulative files processed: %d\n', filesCumulative);
fprintf(fid, 'cumulative hours runtime: %f\n',runtimeCumulative);
fclose(fid);

performance = toc;
avg_per_file = performance / numFiles;
disp("NERD STATS")
stat1 = fprintf("run time: %f s\n",performance);
stat2 = fprintf("runtime average per file: %f s\n",avg_per_file);
disp("LIFETIME NERD STATS")
stat3 = fprintf("total files analyzed: %d\n",filesCumulative);
stat4 = fprintf("total run time: %f hrs\n",runtimeCumulative);
end