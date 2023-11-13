function [data] = optoStimEpoch(data,Hz)
% Separates the optical stimulation TTLs stored in Pe1/ by frequency

if Hz == 10
    var = 100;
elseif Hz == 20
    var = 50;
elseif Hz == 40
    var = 25;
else
    disp('Invalid Hz paramater!')
end


index = (data.epocs.Pe1_.data == var);
data.epocs.Pe1_.onset = data.epocs.Pe1_.onset(index,:);
data.epocs.Pe1_.offset = data.epocs.Pe1_.offset(index,:);



end