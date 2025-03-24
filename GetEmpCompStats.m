

    for k = 1:5
    MDLDATA = [];
    MDLDATA_IND = [];
    for j = 1:10
        data2add = BestFit{j,k}(:,k);
        MDLDATA = [MDLDATA; data2add];
        MDLDATA_IND = [MDLDATA_IND; ones(size(data2add)).*(j)];
    end

    for j = 1:10
    FITrange{k}(j,:) = [min(MDLDATA(MDLDATA_IND==j)) max(MDLDATA(MDLDATA_IND==j))];
    % if k == 1
    % for j = 1:10
    % maxKS_range(j,:) = [min(MDLDATA(MDLDATA_IND==j)) max(MDLDATA(MDLDATA_IND==j))];
    % end
    % elseif ismember(k,[6 7])
    % for j = 1:10
    % FITrange{k}(j,:) = [min(MDLDATA(MDLDATA_IND==j)) max(MDLDATA(MDLDATA_IND==j))];
    % end
    % end

    CV(j,k) = std(MDLDATA(MDLDATA_IND==j))/mean(MDLDATA(MDLDATA_IND==j));

    meanVal(j,k)  = mean(MDLDATA(MDLDATA_IND==j));
    stdVal(j,k)  = std(MDLDATA(MDLDATA_IND==j));
    end
    end

    Emp_maxKSrange = [min(EmpFitType(:,1)) max(EmpFitType(:,1))];

    for j = 1:10
        for k = 1:5
        Emp_range = [min(EmpFitType(:,k)) max(EmpFitType(:,k))];
        data = BestFit{j,k}(:,k);
        INRange = (data>=Emp_range(1) & data<=Emp_range(2));
        meanINRange(j,k) = mean(INRange);
        end
    end

    MAP_names = {'Spatial','Gene coexpression','Receptor similarity','Laminar similarity','Metabolic connectivity','Haemodynamic connectivity','Electrophysiological connectivity','Temporal similarity','Random similarity','Matching'};

    Emp_maxKSrange = [min(EmpFitType(:,6)) max(EmpFitType(:,6))];

    for j = 1:10
        for k = 1:5
        data = BestFit{j,k}(:,6);
        INRange = (data>=Emp_maxKSrange(1) & data<=Emp_maxKSrange(2));
        meanINRange_edge(j,k) = mean(INRange);
        end
    end

    Emp_maxKSrange = [min(EmpDegCorr) max(EmpDegCorr)];

    for j = 1:10
        for k = 1:5
        data = BestFit{j,k}(:,7);
        INRange = (data>=Emp_maxKSrange(1) & data<=Emp_maxKSrange(2));
        meanINRange_deg(j,k) = mean(INRange);
        end
    end

    T = table(MAP_names',meanINRange);


    meanVal