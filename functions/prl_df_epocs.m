function [data] = prl_df_epocs(data,TTLs)

if TTLs == 1
    
    cueTSA = data.epocs.St1_.onset;
    %Makes separate epocs for different trial outcomes
    if isfield(data.epocs, 'CL1_') == 0
        correct_ts_A = 0;
    else
        correct_ts_A = data.epocs.CL1_.onset;
    end

    if isfield(data.epocs, 'IL1_') == 0
        incorrect_ts_A = 0;
    else
        incorrect_ts_A = data.epocs.IL1_.onset;
    end
    pellet_ts_A = data.epocs.Pe1_.onset;
    correct_rewardedA = zeros(double(height(pellet_ts_A)));
    correct_norewardA = zeros(double(height(pellet_ts_A)));
    incorrect_rewardedA = zeros(double(height(pellet_ts_A)));
    incorrect_norewardA = zeros(double(height(pellet_ts_A)));
    for i = 1:height(correct_ts_A)
        var1 = correct_ts_A(i,:);
        var2 = pellet_ts_A;
        [f, ~] = ismember(var1, pellet_ts_A);
        if f == 1
            correct_rewardedA(i,:) = var1;
        elseif f < 1
            correct_norewardA(i,:) = var1;
        end  
    end
    for i = 1:height(incorrect_ts_A)
        var1 = incorrect_ts_A(i,:);
        var2 = pellet_ts_A;
        [f, ~] = ismember(var1, pellet_ts_A);
        if f == 1
            incorrect_rewardedA(i,:) = var1;
        elseif f < 1
            incorrect_norewardA(i,:) = var1;
        end
    end
    
    correct_rewardedA = nonzeros(correct_rewardedA(:,1));
    correct_norewardA = nonzeros(correct_norewardA(:,1));
    incorrect_rewardedA = nonzeros(incorrect_rewardedA(:,1));
    incorrect_norewardA = nonzeros(incorrect_norewardA(:,1));
    
    if isempty(correct_rewardedA) == 1
        correct_rewardedA = 0;
    end
    if isempty(correct_norewardA) == 1
        correct_norewardA = 0;
    end
    if isempty(incorrect_rewardedA) == 1
        incorrect_rewardedA = 0;
    end
    if isempty(incorrect_norewardA) == 1
        incorrect_norewardA = 0;
    end
    data.epocs.cRewA.name = 'cRewA';
    data.epocs.cRewA.onset = correct_rewardedA;
    data.epocs.cRewA.offset = correct_rewardedA+1;
    data.epocs.cRewA.data = ones(height(correct_rewardedA));
    data.epocs.cNoRewA.name = 'cNoRewA';
    data.epocs.cNoRewA.onset = correct_norewardA;
    data.epocs.cNoRewA.offset = correct_norewardA+1;
    data.epocs.cNoRewA.data = ones(height(correct_norewardA))*2;
    data.epocs.iRewA.name = 'iRewA';
    data.epocs.iRewA.onset = incorrect_rewardedA;
    data.epocs.iRewA.offset = incorrect_rewardedA+1;
    data.epocs.iRewA.data = ones(height(incorrect_rewardedA))*3;
    data.epocs.iNoRewA.name = 'iNoRewA';
    data.epocs.iNoRewA.onset = incorrect_norewardA;
    data.epocs.iNoRewA.offset = incorrect_norewardA+1;
    data.epocs.iNoRewA.data = ones(height(incorrect_norewardA))*4;

elseif TTLs == 2
    cueTSC = data.epocs.St2_.onset;
    %Makes separate epocs for different trial outcomes
    if isfield(data.epocs, 'CL2_') == 0
        correct_ts_C = 0;
    else
        correct_ts_C = data.epocs.CL2_.onset;
    end
    
    if isfield(data.epocs, 'IL2_') == 0
        incorrect_ts_C = 0;
    else
        incorrect_ts_C = data.epocs.IL2_.onset;
    end
    pellet_ts_C = data.epocs.Pe2_.onset;
    correct_rewardedC = zeros(double(height(pellet_ts_C)));
    correct_norewardC = zeros(double(height(pellet_ts_C)));
    incorrect_rewardedC = zeros(double(height(pellet_ts_C)));
    incorrect_norewardC = zeros(double(height(pellet_ts_C)));
    for i = 1:height(correct_ts_C)
        var3 = correct_ts_C(i,:);
        var4 = pellet_ts_C;
        [g, ~] = ismember(var3, pellet_ts_C);
        if g == 1
            correct_rewardedC(i,:) = var3;
        elseif g < 1
            correct_norewardC(i,:) = var3;
        end  
    end
    for i = 1:height(incorrect_ts_C)
        var3 = incorrect_ts_C(i,:);
        var4 = pellet_ts_C;
        [g, ~] = ismember(var3, pellet_ts_C);
        if g == 1
            incorrect_rewardedC(i,:) = var3;
        elseif g < 1
            incorrect_norewardC(i,:) = var3;
        end
    end
    correct_rewardedC = nonzeros(correct_rewardedC(:,1));
    correct_norewardC = nonzeros(correct_norewardC(:,1));
    incorrect_rewardedC = nonzeros(incorrect_rewardedC(:,1));
    incorrect_norewardC = nonzeros(incorrect_norewardC(:,1));
    if isempty(correct_rewardedC) == 1
        correct_rewardedC = 0;
    end
    if isempty(correct_norewardC) == 1
        correct_norewardC = 0;
    end
    if isempty(incorrect_rewardedC) == 1
        incorrect_rewardedC = 0;
    end
    if isempty(incorrect_norewardC) == 1
        incorrect_norewardC = 0;
    end
    data.epocs.cRewC.name = 'cRewC';
    data.epocs.cRewC.onset = correct_rewardedC;
    data.epocs.cRewC.offset = correct_rewardedC+1;
    data.epocs.cRewC.data = ones(height(correct_rewardedC));
    data.epocs.cNoRewC.name = 'cNoRewC';
    data.epocs.cNoRewC.onset = correct_norewardC;
    data.epocs.cNoRewC.offset = correct_norewardC+1;
    data.epocs.cNoRewC.data = ones(height(correct_norewardC))*2;
    data.epocs.iRewC.name = 'iRewC';
    data.epocs.iRewC.onset = incorrect_rewardedC;
    data.epocs.iRewC.offset = incorrect_rewardedC+1;
    data.epocs.iRewC.data = ones(height(incorrect_rewardedC))*3;
    data.epocs.iNoRewC.name = 'iNoRewC';
    data.epocs.iNoRewC.onset = incorrect_norewardC;
    data.epocs.iNoRewC.offset = incorrect_norewardC+1;
    data.epocs.iNoRewC.data = ones(height(incorrect_norewardC))*4;
end
