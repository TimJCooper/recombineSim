% Author(s): T.J.Cooper
% Updated: 1/12/2016
% Recombination Event Distributions (NCO/CO/Total)
% Input: Event-summaries, detailing counts of each recombination event type (NCO, CO, NA, Total) per chromosome, and global IED distributions derived from NGS-based post-meiotic analysis of recombination
% Output: Simulated IED distributions (per-cell + population average) under varying conditions including randomised deposition of events or interference (Gamma-derived Hazard-functions) for each event-type.

function RecombineSim(eventfile,varfile,folder,genotype,threshold,samples,mode,COratio,express,alpha,beta)
warning off;

%% Data Import & Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid1 = fopen('ChrSizesS288cH4L2_L2HG.txt','r');
chrsizes = round(cell2mat(textscan(fid1, '%d','HeaderLines',1))/100);
data = readtable(eventfile,'Delimiter','\t');
vars = readtable(varfile,'Delimiter','\t');
vartable = table2array(vars(:,2:4));
indices = find(data{:,{'threshold'}}==threshold & strcmp(data{:,{'Genotype'}},genotype) & ~ismember(data{:,{'type'}}, {'8:0','0:8','0:8_8:0','8:0_0:8'}));
eventl = table2array(data(indices,{'len_mid'}));
data = data(indices,:);
data = sortrows(data,{'Meiosis','chr','midpoint'},{'ascend','ascend','ascend'});
ID = table2array(data(:,{'GenotypeID'}));
uID = unique(ID);
chrindex = table2array(data(:,{'chr'}));
[bincounts,~] = histc(eventl,1:50:round((max(eventl)/50)*1.5)*50);
eventsizes = transpose(1:50*length(bincounts)); resec = 0;
eventw = repelem(bincounts,50);
CO_NCO = (data{:,'CO_NCO'});
eventlist = char(regexprep(CO_NCO, {'NCO','CO','U'},{'A','B','C'}));
eventcount = reshape(crosstab(chrindex(:),eventlist(:)-64,ID(:)),length(chrsizes),[],1);
fclose('all');
if ismember('C',eventlist)==0
    pad = zeros(16,length(uID));
    eventcount = InsertRows(eventcount',pad',[2:2:length(uID)*2])';
end
progress = strcat('Simulating:',genotype,'||','Mode:',mode,'||',folder);
disp(progress);
g = 0;

%% Interference Function & Event Distributions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s=0; w=0;
for m = 1:3:size(eventcount,2)
    s=s+1;
    count = eventcount(:,m:(m-1)+3);
    count(:,3) = sum(count(:,1:3),2);
    if strcmp(express,'Y') == 1
        count(:,1) = 2;
        count(:,3) = 2;
    end
    [~,col] = find(count >= 1);
    IEDnum = sum(count,1)-(histc(col, 1:size(count,2))');
    rowID = find(data{:,{'GenotypeID'}}==uID(s));
    midpoint = table2array(data(rowID,{'midpoint'}))';
    IED = zeros(length(midpoint)-length(unique(chrindex(rowID,:))),3);
    [a,~,subs] = unique([eventlist(rowID,1)-64 chrindex(rowID,1)],'rows');
    [~, I] = sort(subs);
    pos = midpoint(I);
    subs = subs(I,:);
    temp = accumarray(subs,1:numel(subs),[],@(x){abs(diff(pos(x(end:-1:1))))});
    for ii = 1:max(a(:,1))
        vals = [temp{ a(:,1) == ii }];
        IED(1:numel(vals),ii) = vals;
    end
    IED(:,3) = cell2mat(accumarray(chrindex(rowID,:),midpoint(:),[], @(x) {(diff(x))}));
    exp{s} = IED;
    for h=1:3
            w=w+1;
            rIED = IED(:,h);
            rIED(rIED==0)= [];
            gam = fitdist(rIED,'gamma');
            gamfit(w,1) = gam.a; gamfit(w,2) = gam.b;
            ci = paramci(gam,0.05);
            gamfit(w,3:6) = ci(:)';
            int = [];
        if strcmp(mode,'Hazard') == 1
            rIED = IED(:,h);
            rIED(rIED==0)= [];
            if gam.a <= 1.01
                gam.a = 1.02;
            end
            adjgam = gam.b/100;
            haz = Hazard(1:20000,gam.a,adjgam);
            indices = find(haz(1,:)>(1/adjgam)*0.95);
            haz = haz(1,1:min(indices));
            haz = normalize_var(haz,0,1);
            if length(haz) > 5000
                haz = haz(1,1:5000);
            end
            int(1:length(haz)) = fliplr(haz);
            int(length(haz)+1:length(haz)*2) = haz;
        elseif strcmp(mode,'UniHazard') == 1 || strcmp(mode,'MixModel') == 1
            adjgam = beta(h)/100;
            haz = Hazard(1:20000,alpha(h),adjgam);
            indices = find(haz(1,:)>(1/adjgam)*0.95);
            haz = haz(1,1:min(indices));
            haz = normalize_var(haz,0,1);
            if length(haz) > 5000
                haz = haz(1,1:5000);
            end
            int(1:length(haz)) = fliplr(haz);
            int(length(haz)+1:length(haz)*2) = haz;
        elseif strcmp(mode,'Random') == 1
            int = ones(1,10);
        else
            error('Invalid event distribution mode selected. Options: Random, Hazard, UniHazard or MixModel');
        end
        width = length(int)/2;
        cmcount = 0; i = 0;
        intwindows{s,h} = int;
        while cmcount~=samples
            i = i+1;
            dist = cell(1,16);
            rawdist = cell(1,16);
            for j=1:16
                varidx = find(vartable(:,1)==j);
                vardata = vartable(varidx,:);
                varm = zeros(chrsizes(j)*100,1);
                varbuffer = zeros(1000,1);
                varmap = [varbuffer;varm;varbuffer];
                varmap(vardata(:,2)+1000) = 1;
                
                num = count(j,h);
                if num==0
                    continue
                end
                varID = find(vars{:,{'chrom'}}==j);
                coords = table2array(vars(varID,{'pos_c'}))';
                telomereL = round(min(coords)/100);
                telomereR = round(max(coords)/100);
                bound = []; lbound = []; rbound = [];
                model = ones(1,chrsizes(j)-telomereL-(chrsizes(j)-telomereR));
                smodel = ones(1,length(model)*10);
                edgeL = length(model)*5;
                edgeR = length(model)*6;
                for k=1:500
                    pos = edgeL:edgeR;
                    weight = smodel(edgeL:edgeR);
                    classrnd = rand(1,1);
                    if classrnd>(COratio/100)
                        ds = pos(sum(bsxfun(@ge,rand(1,1),cumsum(weight./sum(weight))),2)+1);
                        test = ((ds-(length(model)*5))*100)-50;
                        test = test+1000;
                        vardensity = sum(varmap(test-500:test+500));
                        g = g+1;
                        snpres(g) = vardensity;
                        if vardensity>1000
                            ds = randi([edgeL,edgeR]);
                            bound(k,1:2) = [ds-resec,ds+resec];
                        else
                            smodel(ds-width:ds+width-1)=smodel(ds-width:ds+width-1).*int;
                            bound(k,1:2) = [ds-resec,ds+resec];
                        end
                    else
                        ds = randi([edgeL,edgeR]);
                        bound(k,1:2) = [ds-resec,ds+resec];
                    end
                    bound = sort(bound);
                    lbound = bound(:,1); rbound = bound(:,2);
                    matches = diff([rbound(1:end-1,:) lbound(2:end,:)],[],2)>(threshold/100);
                    a = k-(sum(matches(:)==0));
                    if a==num
                        break
                    else
                    end
                end
                stop = [matches;1];
                start = [1;stop(1:end-1)];
                merge = floor(mean([lbound(start~=0) rbound(stop~=0)],2));
                dist{:,j} = transpose(diff(merge));
                rawIED = transpose(diff(lbound+((rbound-lbound)/2)));
                rawdist{:,j} = rawIED(rawIED>0);
            end
            if sum(cellfun('length',dist))==IEDnum(h)
                cmcount = cmcount+1;
            end
            distdf{s,i,h} = [dist{:}];
            rawdistdf{s,i,h} = [rawdist{:}];
        end
        if i>samples
            idx = find(cellfun('length',distdf(:,:,h))~=IEDnum(h));
            for b=1:length(idx)
                distdf{:,idx(b),h} = [];
                rawdistdf{:,idx(b),h} = [];
            end
        end
    end
end

%% Directories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mkdir(strcat(pwd,'/',folder,'/','Simulations','/',genotype,'/',mode))
chk = exist(strcat(pwd,'/','Results','/','Event_Counts'),'dir');
chk2 = exist(strcat(pwd,'/','Results','/','Experimental'),'dir');
if chk ~= 7 && chk2 ~= 7
    mkdir(strcat(pwd,'/',folder,'/','Event_Counts'))
    mkdir(strcat(pwd,'/',folder,'/','Experimental'))
end

%% Population Output Labels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
types = {'NCO' 'CO' 'Total'};
sets = {1:numel(uID),1:numel(types)};
[p1,p2] = ndgrid(sets{:});
comb = sortrows([p1(:) p2(:)],2);
for e = 1:length(comb)
    AvgPop{e} = strcat(genotype,int2str(uID(comb(e,1))),types{comb(e,2)});
end

%% Results - Population Averaging (Merged)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
resex = permute(distdf, [2 1 3]);
resex = squeeze(mat2cell(resex, size(distdf,2), ones(1, s), ones(1, h)));
stdev = cellfun(@(x) std(sort(cell2mat(x).'), 0, 2), resex, 'un', 0);
resex = cellfun(@(x) [x{:}], resex, 'uniformoutput',false);
resex = cellfun(@sort,resex,'uniformoutput',false);
dec = cellfun(@(x) decimate(x,samples),resex,'uniformoutput',false);
Lmax = max(max(cell2mat(cellfun(@numel,dec,'un',0))));
stdev = cellfun(@(x) [x; nan(max(Lmax(:)) - numel(x), 1)], stdev, 'un', 0);
stdev = cell2mat((stdev(:).'));
b = cellfun(@(c)[c(:);NaN(Lmax-numel(c),1)],dec,'uniformoutput',0);
stre = cell2mat((b(:).'))*100;
T = array2table(stre,'VariableNames',AvgPop);
filename = strcat(pwd,'/',folder,'/','Simulations','/',genotype,'/',mode,'/',genotype,'-',mode,'-',num2str(threshold),'bp','-PopAvg','.txt');
writetable(T,'temp.txt','Delimiter','\t');
replaceinfile('NaN',' ','temp.txt',filename);
delete('temp.txt');

%% Results - Population Averaging (Raw)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rawresex = permute(rawdistdf, [2 1 3]);
rawresex = squeeze(mat2cell(rawresex, size(rawdistdf,2), ones(1, s), ones(1, h)));
rawresex = cellfun(@(x) [x{:}], rawresex, 'uniformoutput', false);
rawresex = cellfun(@sort,rawresex,'uniformoutput',false);
rawdec = cellfun(@(x) decimate(x,samples),rawresex,'uniformoutput',false);
rawLmax = max(max(cell2mat(cellfun(@numel,rawdec,'un',0))));
rawb = cellfun(@(c)[c(:);NaN(rawLmax-numel(c),1)],rawdec,'uniformoutput',0);
rawstre = cell2mat((rawb(:).'))*100;
Tr = array2table(rawstre,'VariableNames',AvgPop);
filename = strcat(pwd,'/',folder,'/','Simulations','/',genotype,'/',mode,'/',genotype,'-',mode,'-',num2str(threshold),'bp','-RawPopAvg','.txt');
writetable(Tr,'temp.txt','Delimiter','\t');
replaceinfile('NaN',' ','temp.txt',filename);
delete('temp.txt');

%% Experimental Output Labels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for e=1:length(uID)
    for r=1:length(types)
        ExpDat{(e-1)*length(types)+r} = strcat(genotype,int2str(uID(e)),types{r});
    end
end

%% Results - Experimental Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = strcat(pwd,'/',folder,'/','Experimental','/',genotype,'-',num2str(threshold),'bp','-ExpIED','.txt');
chk = exist(filename,'file');
if chk ~= 2
    Lmax = max(cellfun('size',exp,1));
    exp = cellfun(@(c) [c;NaN(Lmax-size(c,1),3)],exp,'uniformoutput',0);
    exp = horzcat(exp{:});
    exp(exp==0)=NaN;
    exp = sort(exp);
    T = array2table(exp,'VariableNames',ExpDat);
    writetable(T,'temp.txt','Delimiter','\t');
    replaceinfile('NaN',' ','temp.txt',filename);
    delete('temp.txt');
end

%% Results - Experimental Data (Aggregate)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if chk ~= 2
    filename = strcat(pwd,'/',folder,'/','Experimental','/',genotype,'-',num2str(threshold),'bp','-AggregateExpIED','.txt');
    agg = reshape(permute(reshape(exp,size(exp,1),3,[]),[1,3,2]),[],3);
    agg = sort(agg);
    for y=1:size(agg,2)
        w=w+1;
        idx = ~isnan(agg(:,y));
        gam = fitdist(agg(idx,y),'gamma');
        gamfit(w,1) = gam.a; gamfit(w,2) = gam.b;
        ci = paramci(gam,0.05);
        gamfit(w,3:6) = ci(:)';
    end
    TAgg = array2table(agg,'VariableNames',{'NCOAvg' 'COAvg' 'TotalAvg'});
    writetable(TAgg,'temp.txt','Delimiter','\t');
    replaceinfile('NaN',' ','temp.txt',filename);
    delete('temp.txt');
end

%% Results - Gamma-Fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = strcat(pwd,'/',folder,'/','Experimental','/',genotype,'-',num2str(threshold),'bp','-GammaFit','.txt');
chk = exist(filename,'file');
if chk ~= 2
    GamDat = horzcat(ExpDat,'NCOAggregate','COAggregate','TotalAggregate');
    fit = array2table(gamfit,'RowNames',GamDat,'VariableNames',{'Alpha','Beta','LowerLimA','UpperLimA','LowerLimB','UpperLimB'});
    writetable(fit,'temp.txt','Delimiter','\t','WriteRowNames',true);
    replaceinfile('Row',' ','temp.txt',filename);
    delete('temp.txt');
end

%% Results - Genotype Aggregation (Merged)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
agg = sort(reshape(stre, size(stre,1)*(numel(resex)/3), 3));
filename = strcat(pwd,'/',folder,'/','Simulations','/',genotype,'/',mode,'/',genotype,'-',mode,'-',num2str(threshold),'bp','-AggregateSim','.txt');
TAgg = array2table(agg,'VariableNames',{'NCOAvg' 'COAvg' 'TotalAvg'});
writetable(TAgg,'temp.txt','Delimiter','\t');
replaceinfile('NaN',' ','temp.txt',filename);
delete('temp.txt');

%% Results - Genotype Aggregation (Raw)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
agg = sort(reshape(rawstre, size(rawstre,1)*(numel(rawresex)/3), 3));
filename = strcat(pwd,'/',folder,'/','Simulations','/',genotype,'/',mode,'/',genotype,'-',mode,'-',num2str(threshold),'bp','-RawAggregateSim','.txt');
TAgg = array2table(agg,'VariableNames',{'NCOAvg' 'COAvg' 'TotalAvg'});
writetable(TAgg,'temp.txt','Delimiter','\t');
replaceinfile('NaN',' ','temp.txt',filename);
delete('temp.txt');

%% Event-count Output Labels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
types = {'NCO' 'CO' 'NA' 'Total'};
for e=1:length(uID)
    for r=1:length(types)
        ExpDat{(e-1)*length(types)+r} = strcat(genotype,int2str(uID(e)),types{r});
    end
end

%% Results - Event-count
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = strcat(pwd,'/',folder,'/','Event_Counts','/',genotype,'-',num2str(threshold),'bp','-EventCount','.txt');
chk = exist(filename,'file');
if chk ~= 2
    ecsum = [];
    for q = 1:3:size(eventcount,2)
        ectmp = eventcount(:,q:(q-1)+3);
        ectmp(:,4) = sum(ectmp(:,1:3),2);
        ecsum = horzcat(ecsum,ectmp);
    end
    NCOavg = (round((mean(ecsum(:,1:4:size(ecsum,2)),2))*100)/100);
    NCOstd = (round((std(ecsum(:,1:4:size(ecsum,2)),0,2))*100)/100);
    COavg = (round((mean(ecsum(:,2:4:size(ecsum,2)),2))*100)/100);
    COstd = (round((std(ecsum(:,2:4:size(ecsum,2)),0,2))*100)/100);
    totavg = (round((mean(ecsum(:,4:4:size(ecsum,2)),2))*100)/100);
    totstd = (round((std(ecsum(:,4:4:size(ecsum,2)),0,2))*100)/100);
    ecsum = horzcat(ecsum,NCOavg,NCOstd,COavg,COstd,totavg,totstd);
    sum(ecsum(:,size(ecsum,2)-1));
    for u=1:16
        chrlabels{u,:} = strcat('chr',num2str(u));
    end
    ExpDat = horzcat(ExpDat,'NCOAverage','NCO_STD','COAverage','CO_STD','TotAverage','Tot_STD');
    EC = array2table(ecsum,'RowNames',chrlabels,'VariableNames',ExpDat);
    writetable(EC,'temp.txt','Delimiter','\t','WriteRowNames',true);
    replaceinfile('Row',' ','temp.txt',filename);
    delete('temp.txt');
end

%% Results - Store Database
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = strcat(pwd,'/',folder,'/','Simulations','/',genotype,'/',mode,'/',genotype,'-',mode,'-',num2str(threshold),'bp','.mat');
save(filename)
