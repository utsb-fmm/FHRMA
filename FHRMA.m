function FHRMA(file)
if nargin==0
    [fname,folder]=uigetfile('Examples/*.fhr;*.fhrm;*.dat');
else
    [folder,fname,ext]=fileparts(file);
    fname=[fname ext];
end
load('./FSdataset/expertAnnotations.mat','DB');
n=find(strcmp(fname,{DB.filename}));
if isempty(n)
    s2=inputdlg(sprintf(['When does the second stage start (in min) ?\n' ...
            'Leave empty if no second stage.\n' ...
            'Use 0 if the recording is only second stage.\n' ...
            'For CTG-UHB database, the second stage start at 60 min for all recs.']));
    Stage2Start=str2double(s2{1})*240+1;
else
    Stage2Start=DB(n).Stage2_Start;
end

[FHRDop,FHRScalp,MHR,TOCO]=fhropen([folder fname]);
%MHR=[zeros(1,50) MHR(1:end-50)];
isStage2=zeros(size(FHRDop));
if ~isempty(Stage2Start) && ~isnan(Stage2Start) && isnumeric(Stage2Start)
    isStage2(Stage2Start:end)=1;
end

if any(FHRDop>0) || any(MHR>0)
    [~,~,DopTrue,MHRTrue]=FalseSigDetectDopMHR(FHRDop,MHR,isStage2);
else
    %PDop=zeros(size(FHRDop));
    %PMHR=zeros(size(FHRDop));
    DopTrue=zeros(size(FHRDop));
    MHRTrue=zeros(size(FHRDop));
end
if any(FHRScalp>0)
    [~,ScalpTrue]=FalseSigDetectScalp(FHRScalp,isStage2);
else
    %PScalp=zeros(size(FHRDop));
    ScalpTrue=zeros(size(FHRDop));
end

FHRfull=ScalpTrue;
FHRfull(FHRfull==0)=DopTrue(FHRfull==0);
FHRi=interpolFHR(FHRfull);

[baseline,acc,dec]=aamwmfb(FHRi);
acc=round(acc*4);dec=round(dec*4);

MASigs=[baseline;baseline;baseline];
for j=1:size(acc,2)
    if sum(FHRfull(acc(1,j)+1:acc(2,j))>0)>=10*4 %To validate an A/D we ask tha there is at least 10s of non missing signal
        MASigs(3,acc(1,j)+1:acc(2,j))=FHRi(acc(1,j)+1:acc(2,j));
    end
end
for j=1:size(dec,2)
    if sum(FHRfull(dec(1,j)+1:dec(2,j))>0)>=10*4
        MASigs(2,dec(1,j)+1:dec(2,j))=FHRi(dec(1,j)+1:dec(2,j));
    end
end
MASigs(2,FHRfull==0)=baseline(FHRfull==0);
MASigs(3,FHRfull==0)=baseline(FHRfull==0);
[ TOCOBL,ctu,~ ] = aamwmfb(TOCO*2);TOCOBL=TOCOBL/2;
ctu=round(ctu*4);
CTUSigs=[TOCO;TOCO];
for j=1:size(ctu,2)
    CTUSigs(2,ctu(1,j)+1:ctu(2,j))=0;%TOCOBL(ctu(1,j)+1:ctu(2,j));
end


fhrplot({DopTrue,ScalpTrue,MHRTrue,MASigs,FHRi,FHRDop,FHRScalp,MHR},CTUSigs,['FHRMA - ' fname])

end

