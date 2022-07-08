function EvalFSAllDatasets()
W_MHR=false;WO_MHR=true;
AllStages=0;Stage1=1;Stage2=2;Antepartum=-1;

params={
    'FSDop',AllStages,'DopMHRTrain',W_MHR;
    'FSMHR',AllStages,'DopMHRTrain',W_MHR;    
    'FSScalp',AllStages,'ScalpTrain',W_MHR;
    'FSDop',AllStages,'DopMHRVal',W_MHR;%4
    'FSDop',AllStages,'DopMHRVal',WO_MHR;
    'FSMHR',AllStages,'DopMHRVal',W_MHR;    
    'FSScalp',AllStages,'ScalpVal',W_MHR;
    'FSDop',Antepartum,'DopMHRVal',W_MHR;
    'FSDop',Antepartum,'DopMHRVal',WO_MHR;   
    'FSDop',Stage1,'DopMHRVal',W_MHR;
    'FSDop',Stage1,'DopMHRVal',WO_MHR;
    'FSMHR',Stage1,'DopMHRVal',W_MHR;    
    'FSScalp',Stage1,'ScalpVal',W_MHR;
    'FSDop',Stage2,'DopMHRVal',W_MHR;
    'FSDop',Stage2,'DopMHRVal',WO_MHR;
    'FSMHR',Stage2,'DopMHRVal',W_MHR;    
    'FSScalp',Stage2,'ScalpVal',W_MHR; 
%     'Dop',AllStages,'DopMHRTestCP',W_MHR;%18
%     'Dop',AllStages,'DopMHRTestCP',WO_MHR;
%     'MHR',AllStages,'DopMHRTestCP',W_MHR;   
%     'Dop',Stage1,'DopMHRTestCP',W_MHR;
%     'Dop',Stage1,'DopMHRTestCP',WO_MHR;
%     'MHR',Stage1,'DopMHRTestCP',W_MHR; 
%     'Dop',Stage2,'DopMHRTestCP',W_MHR;
%     'Dop',Stage2,'DopMHRTestCP',WO_MHR;
%     'MHR',Stage2,'DopMHRTestCP',W_MHR;
%     'Scalp',AllStages,'ScalpTest',W_MHR;%27
%     'Scalp',Stage1,'ScalpTest',W_MHR;
%     'Scalp',Stage2,'ScalpTest',W_MHR;
%     'a',AllStages,'DopMHRTestDbS',W_MHR;%30
%     'b',AllStages,'DopMHRTestDbS',W_MHR;
%     'Dop',AllStages,'DopMHRTestDbS',W_MHR%32
    };

 results=[];
 stats=[];

for i=1:length(params)
    disp(params(i,:))
    [s,r]=EvalFSForDataset(params{i,1},params{i,2},params{i,3},params{i,4});
    stats=[stats s]; %#ok<AGROW> 
    results=[results r]; %#ok<AGROW> 
end

save allResutlsFS stats results

f=uifigure;uitable(f,'Data',struct2table(stats),'Unit','Normalized','Position',[0 0 1 1]);