FitTypes = {'maxKS','DegCorr'};

for f = 1:length(FitTypes)
    FitType = FitTypes{f};

    for MdlTypes = 1:3 
    
            if MdlTypes == 1
                    Mdl_names = {'Spatial','Gene coexpression','Receptor similarity','Laminar similarity','Metabolic connectivity','Haemodynamic connectivity','Electrophysiological connectivity','Temporal similarity','Random similarity','Matching'};
                    mdls = 0:9;
                    MdlForms = {'Add','Mult'};
                    Exponents = {'powerlaw','exponential'};
                    save_name_prefix = 'BestMdls_GNM';
                elseif MdlTypes == 2
                    mdls = 1:12;
                    MdlForms = {'Add','Mult'};
                    Exponents = {'powerlaw','exponential'};
                    Mdl_names = cell(1,length(mdls));
                    save_name_prefix = 'BestMdls_TopoMdl_GNM';
                elseif MdlTypes == 3
                    Mdl_names = {'Spatial','Gene coexpression','Receptor similarity','Laminar similarity','Metabolic connectivity','Haemodynamic connectivity','Electrophysiological connectivity','Temporal similarity','Random similarity','Matching'};
                    mdls = 0:9;
                    MdlForms = {'Add'};
                    Exponents = {'exponential'};
                    save_name_prefix = 'BestMdls_WB_GNM';
            end
    
            if ismember(MdlTypes,[1 2])
                mdldata = load('Hansen_networks.mat');
                A_dist = mdldata.A_dist;
                A = mdldata.adj{1};
                d = triu2vec(mdldata.A_dist,1);
                avec = triu2vec(A,1);
                dist_thr{2} = d<30;
                dist_thr{3} = d>=30 & d<=90;
                dist_thr{4}= d>90;
                dist_thr{1}= d>=0;
                d_short = d<30;
                d_mid = d>=30 & d<=90;
                d_long = d>90;
            else
                mdldata = load('Hansen_networks_WB.mat');
                A_dist = mdldata.A_dist;
                A = mdldata.adj{1};
                d = triu2vec(mdldata.A_dist,1);
                avec = triu2vec(A,1);
                dist_thr{2} = d<30;
                dist_thr{3} = d>=30 & d<=90;
                dist_thr{4}= d>90;
                dist_thr{1}= d>=0;
                d_short = d<30;
                d_mid = d>=30 & d<=90;
                d_long = d>90;
            end
    
    
        for formInd = 1:length(MdlForms)
            mdlform = MdlForms{formInd};
        
            for expoInd = 1:length(Exponents)
        
                expo = Exponents{expoInd};
      
                NMdls = length(mdls);
                maxKS = zeros(NMdls,100);
                DegCorr = zeros(NMdls,100);   
                R = zeros(NMdls,100,4);
                EdgeFDR = zeros(NMdls,100,4);
                meanProb = zeros(NMdls,length(avec));
                for mdlIND = 1:NMdls
                    mdl = mdls(mdlIND);
                
                    if MdlTypes == 1
                    if mdl == 9
                        Output = load(['GNM_TopoMdl',num2str(2),'_',mdlform,'_',expo,'.mat']);
                    else
                        Output = load(['GNM_Mdl',num2str(mdl),'_',mdlform,'_',expo,'.mat']);
                    end
                    elseif MdlTypes == 2
        
                        Output = load(['GNM_TopoMdl',num2str(mdl),'_',mdlform,'_',expo,'.mat']);
                        mdlname = Output.Input.TopoType;
                        Mdl_names{mdlIND} = mdlname;
        
                    elseif MdlTypes == 3
                        Output = load(['WB_GNM_Mdl',num2str(mdl),'_',mdlform,'_',expo,'.mat']);
                    end
                
                                   
    
                    maxKS(mdlIND,:) = Output.BestFit.(FitType).maxKS;
                    bnets = Output.BestFit.(FitType).b;
                    DegCorr(mdlIND,:) = Output.BestFit.(FitType).DegCorr;
                    meanProb(mdlIND,:) = Output.BestFit.maxKS.meanProb;
                
                    n = length(A);
                    
                    for i = 1:length(bnets)
                        b = bnets{i};
                        B = zeros(n);
                        B(b) = 1;
                        B = B + B';
                        bvec = triu2vec(B,1);
                        %MdlDeg(mdlIND,:) = MdlDeg(mdlIND,:)+sum(B);
                    
                        for j = 1:4
                            thr = dist_thr{j};
                            athr = avec.*thr;
                            bthr = bvec.*thr;        
                    
                            %r0(mdlIND,i,j) = sum((athr==0).*(bthr==0))./sum(athr==0);
                            R(mdlIND,i,j) = sum(athr.*bthr)./sum(athr);
                    
                            EdgeFDR(mdlIND,i,j) = sum(athr==0&bthr==1)./sum(bthr);
                        end
                    
                    end
                
                end
        
            save([save_name_prefix,'_',FitType,'_',mdlform,'_',expo,'.mat'],'R','EdgeFDR','maxKS','DegCorr','Mdl_names','expo','mdlform','meanProb')
    
            end
        
        end
    
    end

end
