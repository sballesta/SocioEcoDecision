%% WRITTEN BY SEBASTIEN BALLESTA (ballesta@unistra.fr) v04092023
%%% Fig 1

clc
close all
clear all
load ('confdat')

% PARAMETERS
limsessmod = 3; %Filtering
limmksess = 10;
nmodisp = 4; % MODULE DISPONIBLES
nsec1 = 25; nsec2 = 1 ;
mkpref = 1 ;

PLOT1 = 1;

%%% ALGO:
modid2=modid;
modid2(modid2==0)= Inf;
modid2(modid2<=limsessmod)= NaN;

mask=sum(isnan(modid2(:,2:end)),2) + sum(isinf(modid2(:,2:end)),2)>nmodisp;
modid2(mask,:) = [];
mask2=sum(~isinf(modid2(:,2:end)),2)==nmodisp; % MODULE DISPONIBLE
modid2(~mask2,:) = [];

nif=sum(sum(isinf(modid2(:,2:end)),1)>0);

if 1
    while nif ~= 4-nmodisp
        % KEEP THE N MODULE THE MOSTE USED
        tmp=sum(isinf(modid2(:,2:end)),1); tmp(tmp==0)=NaN; [~, id]=min(tmp);
        mask=isinf(modid2(:,id+1));
        modid2(mask,:)= [];
        nif=sum(sum(isinf(modid2(:,2:end)),1)>0);
    end
    selectmods=sum(~isinf(modid2(:,2:end)),1)>0; selectmods = find(selectmods==1);

else
    % selectmods = [1 2 3 4];
end

bestdays=modid2(:,1);
disp(['Selected Modules: ' num2str(selectmods)])
disp(['Selected Days: ' num2str(numel(bestdays))])

% FILTERING DATA
dates2=datenum(data.Date); mask=[];
for n=1:numel(bestdays)
    mask(n,:)=dates2==bestdays(n);
end
datemask=logical(sum(mask));
data2=data(datemask,:);

allelo=readtable('elo_matrix_All.xlsx');
prefmod=[]; D=0; elolo=[];
alldate = (unique(data2.Date));
stepdate=round(numel(alldate)/4);

for dte=1 % dte=1:stepdate:numel(alldate)
    D=D+1;
    disp(dte)
    data3=data2;

    %%%% PREFERENCE ALONE FOR ALL MONKEYS
    [rez, MDAY2, TSTART2, MKNAME2, MKMOD2, MKID, MKNAMEID, CDAY2] = dodays_v2(data3,0);

    rezmod=[];
    for k=1:size(rez,1)
        rezmod=[rezmod; rez{k,2}];
    end

    bestmod=[];
    for n=1:4
        bestmod(1,n)=sum(rezmod==n);
    end

    % bestmod(2,:)=bestmod./sum(bestmod)*100;
    theomod=[];
    for n=1:4
        theomod(n)=round(sum(bestmod(1,:)).*1/nmodisp);
    end

    % CHI SQUARE!
    % Chi-square test, by hand
    observed = bestmod;
    expected = theomod;
    chi2stat = sum((observed-expected).^2 ./ expected);
    p = 1 - chi2cdf(chi2stat,1);

    %%% PREFERENCE BY MONKEY
    rezname=mkfinal;
    mkmods=zeros(numel(rezname),4);

    for k=1:size(rez,1)
        for mk=1:size(rezname,1)
            for z=1:size(rez{k,1},1)
                if strcmp(rez{k,1}(z),rezname{mk,1})
                    mkmods(mk,rez{k,2}(z)) = mkmods(mk,rez{k,2}(z)) +1;
                end
            end
        end
    end

    mask=sum(mkmods,2)<limmksess;
    mkmods(mask,:)=NaN;
    pmk=mkmods./sum(mkmods,2)*100;

    if PLOT1

        mkmods2=mkmods(:,selectmods);
        myC= [0 0 .8
            1 0.6 0
            .4 0.4 .4
            0 0.8 1
            ];
        figure(4) %%%% TO DO: ADD STAT CHISQUARE PER MONKEY !!!
        h = bar((mkmods2./sum(mkmods2,2)*100),'stack'); hold on; box off;
        plot([-1 28], [100/numel(selectmods) 100/numel(selectmods)],'k--','HandleVisibility','off');
        for k=1:numel(selectmods)
            set(h(k),'facecolor',myC(k,:))
        end

        ax=legend(h, {'MALT #1', 'MALT #2', 'MALT #3', 'MALT #4'}, 'Location','Best','FontSize',8);
        xticks(1:size(rezname,1)); xticklabels(rezname); xtickangle(45)
        yticks([0:25:100])
        ylabel('% of MALT usage')
        disp( ['Total trials # : ' num2str(nansum(nansum(mkmods))) ] )
        axis tight
        ax = gca;
        ax.FontSize = 16;

    end

    % CALCULATE RELATIVE VALUE FOR EACH MONKEY
    [pmkval, pmkid]=sort(pmk,2,'descend');
    pmkdiff=abs(diff(pmkval,1,2));
    prefmod{D}=mkmods./sum(mkmods,2);

end


%%%%% DO CHI2
pjuice = [];
for mk=1:size(mkmods,1)
    % Observed data
    n1 = mkmods(mk,:); N1 = sum(n1);
    n2 = round([0.25 0.25 0.25 0.25].*sum(n1)); N2 = sum(n1);
    % Pooled estimate of proportion
    p0 = (n1+n2) / (N1+N2)
    % Expected counts under H0 (null hypothesis)
    n10 = N1 * p0;
    n20 = N2 * p0;
    % Chi-square test, by hand
    observed = [n1 N1-n1 n2 N2-n2];
    expected = [n10 N1-n10 n20 N2-n20];
    [h,pjuice(mk),stats] = chi2gof([1:16],'freq',observed,'expected',expected,'ctrs',[1:16],'nparams',2)
end

% save('conflictrez')

%% FIG S1 (or FIG 3a if allmk = 27)
clear all
clc
close all

load ('atab_4.mat')

allmk =  1:26; %% % (27 for ALL MONKEY analysis = fig 3a sinon = fig S1)
% allmk = 27; %%

c=colormap('copper');
stab=table;
nskip=0;
PLOT=1;
% FILTER FOR DIFFERENT SOCIAL SITUATIONS (see code)
condfilter = 0 ; %%



for nmodisp = [ 4 ]
    disp('... building the GLM ')

    mkrank2=mkidx.rank';

    deltarank{nmodisp} = ones( numel(mkrank2) ,6) .*NaN ;
    steepness{nmodisp} = ones( numel(mkrank2) ,6) .*NaN ;
    pref{nmodisp} = ones( numel(mkrank2) ,6) .*NaN ;
    pval2D{nmodisp} = ones( numel(mkrank2) ,6) .*NaN ;


    for actselect =  allmk % (27 for ALL MONKEY analysis = fig 3a otherwise = fig S1)
        try
            disp(['...analysing monkey: ' mkfinal{actselect}])
        catch
            disp(['...analysing ALL monkeys '])
        end

        load ('atab_4.mat','atab')

        % SELECT ONE ACTORS
        if actselect<27
            mask=atab.ActorID~=actselect;
            atab(mask,:)=[];

        elseif actselect==27
        end

        % CLEAN atab of previous SELF PRESENCE or Excluded Actor
        mask=(atab.ActorID==0 | isnan(atab.ActorID));
        atab(mask,:)=[];

        if 1 %%% CALCULATE MODULE PREFERENCE IN THIS DATASETS FOR EACH MONKEY and PATCH

            mask=(atab.OccupMod == 0);
            atab2=atab(mask,:);
            ratiopref=zeros(numel(mkrank2),4).*NaN;
            selfpref=atab2.ChosenMod==selectmods;
            for mk=1:numel(mkrank2)
                mask=atab2.ActorID==mk;
                mkprefer=selfpref(mask,:);
                ratiopref(mk,selectmods)=sum(mkprefer)/sum(sum(mkprefer));
            end

            for z=1:size(atab,1)
                idmk=atab.ActorID(z);
                for mod=1:4
                    eval(['atab.SelfPrefMod' num2str(mod) '(z) = ratiopref(idmk,mod);' ])
                end
            end
        end

        % CLEAN situation when MALT both empty
        mask=(atab.OccupMod ==0);
        emptysess=sum(mask);
        if condfilter >=0
            atab(mask,:)=[];
        end

        if condfilter==1
            % Remove if at least one empty mod
            disp('Remove if at least one empty mod ')
            mask=atab.OccupMod<numel(selectmods);
            atab(mask,:)=[];
        elseif condfilter==2
            % Select only at least one empty mod
            mask=atab.OccupMod==numel(selectmods);
            atab(mask,:)=[];
        elseif condfilter==-1
            mask=(atab.OccupMod >0);
            atab(mask,:)=[];
        end


        % FILTER TABLE FOR PARTICULAR INDIVIDUALS
        allmod = combnk(1:4,2);
        %ANALYSIS OF ATAB
        for k = 1:size(atab,1)
            atab.pairmod(k,:)=zeros(1,6);
            for i=1:size(allmod,1)
                check=[];
                for j=1:numel(selectmods)
                    check(j,:)=atab.DispMod(k,j)==allmod(i,:);
                end

                if sum(sum(check))>=2 && ( atab.ChosenMod(k)==allmod(i,1) || atab.ChosenMod(k)==allmod(i,2) )
                    atab.pairmod(k,i)=1;
                    eval([ 'atab.prefdiff(k,i) =  [ round(atab.SelfPrefMod' num2str(allmod(i,1)) '(k),4) -  round(atab.SelfPrefMod' num2str(allmod(i,2)) '(k),4)] ; '])
                    eval([ 'atab.rankdiff(k,i)= round((atab.RankMod' num2str(allmod(i,1)) '(k)) - (atab.RankMod' num2str(allmod(i,2)) '(k)),2);' ]);
                else
                    atab.prefdiff(k,i)=NaN;
                    atab.rankdiff(k,i)=NaN;
                end
            end
        end

        if ~isempty(atab)
            for i=1:6

                if actselect~=27
                    maskm=atab.pairmod(:,i)==1;
                elseif actselect ==27
                    maskm=~isnan(atab.pairmod(:,i)>0);
                end

                if sum(maskm)==0 || (actselect==27 && i>1)
                    continue
                end
                temptab=atab(maskm,:);
                %%% MODELLING
                try
                    ttype =  [ temptab.prefdiff(:,i) temptab.rankdiff(:,i)  ];
                catch
                    %         disp('error')
                    continue
                end
                trialtypes=unique(ttype,'rows');
                mask=sum(isnan(trialtypes),2)>0; % REMOVE NaN
                trialtypes(mask,:)=[];


                for k = 1:size(temptab,1)
                    tt=false(size(trialtypes));
                    for z=1:size(trialtypes,1)
                        tt(z,:) = tt(z,:) + trialtypes(z,:)==ttype(k,:);
                    end
                    tt=sum(tt,2);
                    tt=find(tt==size(ttype,2));

                    if ~isempty(tt)
                        temptab.Trialtype(k)=tt;
                    else
                        temptab.Trialtype(k)=NaN;
                    end
                end

                mask=isnan(temptab.Trialtype);
                temptab(mask,:)=[];
                choice=[]; binosize=[];
                for n=1:size(trialtypes,1)
                    mask=temptab.Trialtype==n;
                    choice(n) = nanmean(temptab.ChosenMod(mask)==allmod(i,1));
                    choice(n) = nansum(temptab.ChosenMod(mask)==allmod(i,1));
                    binosize (n) = sum(mask);
                end

                gtab=table;
                prefmod=trialtypes(:,1);
                gtab.rankmod=zscore(-trialtypes(:,2));
                gtab.choice=[ choice' binosize' ];
                mask=binosize<=nskip;
                gtab(mask,:)=[];
                prefmod(mask,:)=[];

                if sum(~mask)>1

                    if PLOT==1
                        figure(1);
                        if actselect ~=27
                            subplot(4,7,actselect);
                        end

                        if actselect == 27
                            list_rank=unique (gtab.rankmod);
                            gtab_temp=table;
                            for rk=1:numel(list_rank)
                                idr = list_rank(rk);
                                mask=gtab.rankmod==idr;
                                gtab_temp.rankmod(rk)=idr;
                                gtab_temp.choice(rk,:)=sum(gtab.choice(mask,:));
                            end
                            gtab=gtab_temp;
                            plot(gtab.rankmod,gtab.choice(:,1)./gtab.choice(:,2),'o','MarkerSize', 5,'MarkerEdgeColor', c((i-1)*45+1,:),'HandleVisibility','off');
                        else
                            plot(gtab.rankmod,gtab.choice(:,1)./gtab.choice(:,2),'.','MarkerSize', .5,'MarkerEdgeColor', c((i-1)*45+1,:),'HandleVisibility','off');
                        end

                        ylim([-.1 1.1])
                        axis square
                        if  actselect == 27 || actselect == 1
                            ylabel( 'p of chosing MALT X  vs MALT Y' );
                            xlabel('Rank Monkey X - Rank Monkey Y');
                        end
                        hold on
                    end

                    glmrez = fitglm(gtab,'choice ~ rankmod ','Distribution','Binomial','Link','logit','BinomialSize',binosize(~mask))


                    if PLOT==1
                        if size(gtab,1)>2
                            %%%%%%% FINDING INDIFFERENCE PTS
                            a0=glmrez.Coefficients.Estimate(1);
                            a1=glmrez.Coefficients.Estimate(2);
                            hogit= log(.5/(1-.5));
                            deltarank{nmodisp}(actselect,i) = -a0/a1;
                            steepness{nmodisp}(actselect,i) = a1 ;
                            pref{nmodisp}(actselect,i) = prefmod(1);
                            pval2D{nmodisp}(actselect,i) = glmrez.Coefficients.pValue(2);

                            if glmrez.Coefficients.pValue(2)<=0.05
                                li='-';
                                si=1.5;
                            else
                                li='--';
                                si=.5;
                            end

                            figure(1)
                            if actselect~=27
                                subplot(4,7,actselect)

                                title( mkfinal{actselect})
                                if abs(-a0/a1) <100
                                    plot( -a0/a1 , .5 ,'o','MarkerEdgeColor', c((i-1)*45+1,:),'HandleVisibility','off');
                                    plot( -a0/a1 , .5 ,'+','MarkerEdgeColor', c((i-1)*45+1,:),'HandleVisibility','off');
                                end

                            else
                                title( 'All monkeys')
                            end

                            % MANUAL SELECTION OF X AXIS LIMIT for each
                            % monkey
                            if sum(actselect == [ 11 13 14 16 20 21 22 25]) ==1 ; limt=5;
                            elseif  sum(actselect == [6 17 23 24]) ==1 ;  limt=10;
                            else; limt=2.5;  end

                            ypred=predict( glmrez, (-limt:.1:limt)' );
                            plot((-limt:.1:limt)',ypred,'LineWidth',si,'LineStyle',li,'Color',c((i-1)*45+1,:));
                        else
                            title( mkfinal{actselect})
                        end
                    end
                    %                     plot([-28:1:28]',ypred,'b');

                else

                    if PLOT==1
                        figure(actselect)
                        title( mkfinal{actselect})
                        axis off
                    end

                    deltarank{nmodisp}(actselect,i) = NaN;
                    steepness{nmodisp}(actselect,i) = NaN ;
                    pref{nmodisp}(actselect,i) = NaN;
                    pval2D{nmodisp}(actselect,i) = NaN;
                end

                if actselect ~=27
                    gtab.prefmod=prefmod; %%% ZSCORE ?
                    gtab.actid=ones(size(gtab,1),1)*actselect;
                    gtab.modisp=ones(size(gtab,1),1)*nmodisp;
                    gtab.modid=ones(size(gtab,1),2).*allmod(i,:);
                end

                stab=[ stab; gtab ];

            end
        else
            disp('...empty tab...')
        end
        box off
        drawnow

        if actselect >=26
            legend ('MALT 1 vs 2', 'MALT 1 vs 3', 'MALT 1 vs 4', 'MALT 2 vs 3', 'MALT 2 vs 4', 'MALT 3 vs 4')
        end

    end

end


%% (Figure 4) INDIFFERENCE POINT AND STEEPNESS ANALYSIS
close all
clc

if ~exist("stab");
    clear all
    load mkid
    load mkrank_final
    load stab
end

for nmodisp=4
    allrez=[];
    for actselect=1:numel(mkrank2)%
        % TO UNCOMMENT FOR PJUICE
        %     for actselect=find(pjuice<0.001)%
        for i=1:6
            allrez=[allrez; [deltarank{nmodisp}(actselect,i) steepness{nmodisp}(actselect,i) pref{nmodisp}(actselect,i)  pval2D{nmodisp}(actselect,i) mkrank2(actselect) pjuice(actselect) ] ];
        end
    end
end

%%% FILTER ALLREZ
allrez_out=[];
mask = (allrez(:,4))>.05;
allrez_out = allrez(mask,:);
allrez(mask,:)=[];

allrez_pjuice = [];
maskpjuice = allrez(:,6)<.01
allrez_pjuice = allrez(maskpjuice,:);
allrez(maskpjuice,:)=[];

figure(1)
subplot(1,2,1)
plot(allrez(:,1) , allrez(:,3) ,'ko','Markersize',5.5);  hold on
plot(allrez_pjuice(:,1) , allrez_pjuice(:,3) ,'bx','Markersize',5.5);
xlabel('Rank indifference point'); ylabel('Difference in MALT preferences');
axis([-5 5 -1 1 ]); axis square; box off; lsline;
plot(allrez_out(:,1) , allrez_out(:,3) ,'ks','Markersize',4);
[R, P]=corrcoef(allrez(:,1) , allrez(:,3));
[R2, P2]=corrcoef(allrez_pjuice(:,1) , allrez_pjuice(:,3))
title( { ['R = ' num2str(round(R2(2),3)) ' p = ' num2str(round(P2(2),3)) ] ...
    ['R = ' num2str(round(R(2),3)) ' p = ' num2str(round(P(2),3)) ]} )


subplot(1,2,2)
plot((allrez(:,2)) , abs(allrez(:,3)) ,'ko','Markersize',4); hold on;
plot((allrez_pjuice(:,2)) , abs(allrez_pjuice(:,3)) ,'bx','Markersize',4);
lsline;
xlabel('Steepness'); ylabel('Absolute difference in MALT preferences');
[R, P]=corrcoef((allrez(:,2)) , abs(allrez(:,3)))
[R2, P2]=corrcoef((allrez_pjuice(:,2)) , abs(allrez_pjuice(:,3)))
title( { ['R = ' num2str(round(R2(2),3)) ' p = ' num2str(round(P2(2),3)) ] ...
    ['R = ' num2str(round(R(2),3)) ' p = ' num2str(round(P(2),3))] } )
axis([-0.5 3 0 1 ]); axis square; box off; lsline
plot(allrez_out(:,2) , abs(allrez_out(:,3)) ,'ks','Markersize',3.5);

%
% scatter( nanmean(allrez(:,2)), nanmean(abs(allrez(:,3))),100,'MarkerEdgeColor','w',...
%         'MarkerFaceColor',[.5 .5 .5] );
% scatter( nanmean(allrez_pjuice(:,2)), nanmean(abs(allrez_pjuice(:,3))),100,'MarkerEdgeColor','w',...
%         'MarkerFaceColor',[.1 .1 .8] );


figure(2)
h = bar ([nanmean(allrez(:,2)), nanmean(allrez_pjuice(:,2))] ) ; hold on;
errorbar(h.XData, h.YData, [ nanstd(allrez(:,2))/sqrt(numel(allrez(:,2))) , nanstd(allrez_pjuice(:,2))/sqrt(numel(allrez_pjuice(:,2))) ],'k.' )
[p, H]=ranksum((allrez(:,2)),(allrez_pjuice(:,2)))
ylabel('Steepness')



%% Figure 5a STAB ANALYSIS
close all
clc
warning off

if ~exist("stab")
    clear all
    load mkid
    load mkrank_final
    load stab
end

mkrank2=mkidx.rank';

[~ , idx]=sort(mkrank2);

Z=0;

%          T = [2 2; 1 0; 0 1;]; %FULL MODEL
%            T = [0 0; 1 0; 0 1; 1 1;]; % WITH INTERACTION
T = [0 0; 1 0; 0 1;]; % NO INTERACTION
rezglm=zeros(26,6).*NaN; pf=[]; rf=[];  hie=0;

% 27 is all monkey
mkselection=27;
%  mkselection=1:26;
%  mkselection=find(pjuice<0.01);

for actselect=  mkselection
    if actselect<27
        Z=Z+1;
        mask=stab.actid==actselect;
        tetab=stab(mask,:);
        tetab.prefmod=zscore(tetab.prefmod);
        pf(actselect)=mean(tetab.prefmod);
        rf(actselect)=mean(tetab.rankmod);

    elseif actselect == 27

        tetab=stab;
        tetab.rankmod = round(tetab.rankmod,1);
        tetab.prefmod = round(tetab.prefmod,1);

        % GATHER DATA FOR ALL ACTORS
        ctype=unique ([tetab.rankmod tetab.prefmod ],'rows' );
        poptab=table;
        for i=1:size(ctype,1)
            mask = tetab.rankmod==ctype(i,1) & tetab.prefmod==ctype(i,2);

            if sum(sum(tetab.choice(mask,:),2)) >= 10
                poptab.rankmod(i) = ctype(i,1) ;
                poptab.prefmod(i) = ctype(i,2) ;
                poptab.choice(i,:) = sum(tetab.choice(mask,:),1);
            else
                poptab.rankmod(i) = ctype(i,1) ;
                poptab.prefmod(i) = ctype(i,2) ;
                poptab.choice(i,:) = [NaN NaN];
            end
        end

        % FILTER OUT
        tetab=poptab;
        tetab.prefmod=zscore(tetab.prefmod);
        tetab.rankmod=zscore(tetab.rankmod);

    end

    if actselect<=27
        disp(actselect)
    else
        disp( poolact)
    end

    try

        tetab(isnan(tetab.choice(:,2)),:)=[];
        glmrez = fitglm ( [tetab.rankmod tetab.prefmod ],tetab.choice,T,'VarNames',{'Rank', 'PrefDiff','Choice'},'Distribution','Binomial','Link','logit' )

        if actselect~=27
            rezglm(actselect,:)= [ glmrez.Coefficients.Estimate' glmrez.Coefficients.pValue' ];
        end

    catch
        if actselect~=27
            disp('...nodata')
            try
                rezglm(actselect,:)= [ NaN NaN NaN NaN NaN NaN ];
            catch
                rezglm(actselect,:)= [ NaN NaN NaN NaN NaN NaN NaN NaN ];
            end
        end
    end

    %%%%%% FIGURE
    if PLOT ==1
        %         figure(actselect)
        figure(667)
        %  subplot(4,7,actselect)
        if actselect<27
            subplot(4,7,Z)
        end

        if actselect<27
            try
                title({ [ mkfinal{actselect} '  Beta Rank = ' num2str(round(glmrez.Coefficients.Estimate(2),3)) ' p = ' num2str(round(glmrez.Coefficients.pValue(2),3))]; ...
                    ['Beta Juice = ' num2str(round(glmrez.Coefficients.Estimate(3),3)) ' p = ' num2str(round(glmrez.Coefficients.pValue(3),3))] ; ...
                    ['Beta Juice*Rank = ' num2str(round(glmrez.Coefficients.Estimate(4),3)) ' p = ' num2str(round(glmrez.Coefficients.pValue(4),3))] ;
                    }) ; hold on;
            catch
                title({ [ mkfinal{actselect} 'Beta Rank = ' num2str(round(glmrez.Coefficients.Estimate(2),3)) ' p = ' num2str(round(glmrez.Coefficients.pValue(2),3))]; ...
                    ['Estimate Juice = ' num2str(round(glmrez.Coefficients.Estimate(3),3)) ' p = ' num2str(round(glmrez.Coefficients.pValue(3),3))] ; ...
                    }) ; hold on;
            end

        else
            try
                title({ [' All monkeys :  Beta Rank = ' num2str(round(glmrez.Coefficients.Estimate(2),3)) ' p = ' num2str(round(glmrez.Coefficients.pValue(2),3))]; ...
                    ['Beta Juice = ' num2str(round(glmrez.Coefficients.Estimate(3),3)) ' p = ' num2str(round(glmrez.Coefficients.pValue(3),3))] ; ...
                    ['Beta Juice*Rank = ' num2str(round(glmrez.Coefficients.Estimate(4),3)) ' p = ' num2str(round(glmrez.Coefficients.pValue(4),3))] ;
                    }) ; hold on;
            catch
                title({ [ ' All monkeys : Beta Rank = ' num2str(round(glmrez.Coefficients.Estimate(2),3)) ' p = ' num2str(round(glmrez.Coefficients.pValue(2),3))]; ...
                    ['Estimate Juice = ' num2str(round(glmrez.Coefficients.Estimate(3),3)) ' p = ' num2str(round(glmrez.Coefficients.pValue(3),3))] ; ...
                    }) ; hold on;
            end
        end


        x2 = tetab.prefmod ;
        y2 = tetab.rankmod ;
        z2 = tetab.choice(:,1)./tetab.choice(:,2) ;

        S=ones(size(x2)).*20;
        y=-2.2:.33:2.2;
        x=-2.2:.33:2.2;


        x_grid = reshape(repmat(x, numel(x), 1), 1, numel(x) ^ 2);
        y_grid = repmat(y, 1, numel(y));
        % Plot the data
        col1=unique(z2);
        col2=gray(round(numel(col1)*1.05));

        col=zeros(numel(z2),3);
        for n=1:numel(col1)
            mask=z2==col1(n);
            col(mask,:)=col(mask,:) + col2(n,:);
        end
        pl = scatter3(x2,y2,z2,S,col,'o','filled','MarkerEdgeColor',[.5, .5, .5]);
        z= reshape(predict(glmrez, [x_grid; y_grid]'), numel(x), numel(y))';
        % Plot the model on the scatter plot
        hold on;
        s = mesh(x, y, z,'FaceAlpha',.5,'EdgeColor', [.5, .5, .5], 'EdgeAlpha',1);

        if  actselect >=27
            colorbar
        end

        colormap gray;

        if hie>0
            surfhie{hie,1}=x;
            surfhie{hie,2}=y;
            surfhie{hie,3}=z;
            surfhie{hie,4}=glmrez;

            surfhie{hie,5}=x2;
            surfhie{hie,6}=y2;
            surfhie{hie,7}=z2;
        end

        if 0
            c = colorbar;
            c.Label.String = 'p of chosing Module#1';
        end

        zlabel('p of chosing MALT X');
        ylabel( 'Rank MALT X - MALT Y' );
        xlabel( 'Juice preference MALT X - MALT Y' );
        axis square
        view(45,30)

    end

end

warning on

%%
%Figure 5BC
close all
clear all
clc
warning off
load rezglm
load mkid
load mkrank_final
load stab
mkidy = readtable ('mkid_Socioeco.xlsx');

rg=rezglm;

% SEXE & ID & PJUICE
[ ~, tdx] = sort(mkidy.idname); mkidy=mkidy(tdx,:);
mkidx.maltid = mkidy.id;
mkidx.sexe=mkidy.sexe;

mkg=mkidx.sexe==1;

figure(127)
subplot(1,2,1)
plot(mkidx.rank(mkg),rg(mkg,2),'bs','Markersize',5);  hold on;
plot(mkidx.rank(~mkg),rg(~mkg,2),'rx','Markersize',5);lsline; axis square;
xlabel('Monkey Rank'); ylabel('Beta of Other Rank');
[R P]=corr(rg(mkg,2),mkidx.rank(mkg),'Type','Spearman');
[R2 P2]=corr(rg(~mkg,2),mkidx.rank(~mkg),'Type','Spearman');
[Rall Pall]=corr(rg(:,2),mkidx.rank,'Type','Spearman');
title( { ['R = ' num2str(round(R(),3)) ' p = ' num2str(round(P(),3)) ] ...
    ['R = ' num2str(round(R2(),3)) ' p = ' num2str(round(P2(),3)) ] }) ;
box off

subplot(1,2,2)
%     plot(mkidx.rank,rg(:,3),'ko','Markersize',4);  hold on; lsline; axis square;
plot(mkidx.rank(mkg),rg(mkg,3),'bs','Markersize',5);  hold on;
plot(mkidx.rank(~mkg),rg(~mkg,3),'rx','Markersize',5);lsline; axis square;
xlabel('Monkey Rank'); ylabel('Beta Juice Pref');
[R, P]=corr(mkidx.rank,rg(:,3),'Type','Spearman');
title( ['R = ' num2str(round(R(),3)) ' p = ' num2str(round(P(),3)) ] ) ;
box off

%% FIGURE S5 And stepwise analysis
%     keyboard
clc
close all
clear all

load mkid
load mkrank_final
load stab
load rezglm

mkrank2=mkidx.rank';

mkidy = readtable ('mkid_Socioeco.xlsx');
load ('prefmalt')

pfmod = [];


for j = 1: size(prefmod{1},1)
    tmp=[];
    for i=1:4

        tmp(i)= abs(prefmod{1}(j,i)-0.25);
    end
    pfmod(j) = nanmean(tmp);
end


mkidy.pjuice =pjuice'; % ADDING PJUICE TO MATRIX

% FILTER MKEY WITHOUT PREF
% pmask=logical(sum(mkidy.id==find(pjuice>0.05),2));
% mkidy(mask,:)=NaN;
% mkidx(mask,:)=NaN;
% mkrank2(mask)=NaN;

% MALTID & JUICE


% SEXE & ID & PJUICE
[ ~, tdx] = sort(mkidy.idname); mkidy=mkidy(tdx,:);
mkidx.maltid = mkidy.id;
mkidx.sexe=mkidy.sexe;
mkidx.pjuice=pfmod';

% AGE
mkidx.age=mkidy.age-0.5;
mkelo = mkidx.elo  ;
mksexe=mkidx.sexe;
mkage=mkidx.age;
mkna=mkidy.idname;
mask=~isnan(rezglm(:,6)) & sum(rezglm,2)>0;
rg=rezglm; rg(~mask,:)=rg(~mask,:).*NaN;


% DOSTEPWISE? % CAREFULL WITH THIS effect of gender but questionnable
glmtab=table;
glmtab.age=zscore(mkidx.age);
glmtab.sexe=zscore(mkidx.sexe);
glmtab.rank=zscore(mkidx.rank);
glmtab.pjuice=zscore(mkidx.pjuice);
glmtab.betarank=rg(:,2);
grank=stepwiseglm(glmtab,'betarank ~ rank + sexe + age + pjuice ','Upper','interactions','Criterion','bic')


glmtab.residuals=grank.Residuals.Pearson;
glmtab.id = mkna;
[~,idx]=sort(glmtab.residuals)


glmtab_rank=glmtab(idx,:);
% glmtab_rank.idn=[26:-1:1]';
glmtab_rank.idn=[1:26]'
mks=glmtab_rank.sexe>1;
plot(glmtab_rank.residuals(mks),glmtab_rank.idn(mks),'rx'); hold on
plot(glmtab_rank.residuals(~mks),glmtab_rank.idn(~mks),'bs');
text(glmtab_rank.residuals+ .01, glmtab_rank.idn, glmtab_rank.id)
box off; axis square;

%     fitglm(glmtab,'betarank ~ rank * sexe ')

glmtab=table;
glmtab.age=zscore(mkidx.age);
glmtab.sexe=zscore(mkidx.sexe);
glmtab.rank=zscore(mkidx.rank);
glmtab.pjuice=zscore(mkidx.pjuice);
glmtab.betajuice=rg(:,3);
gjuice=stepwiseglm(glmtab,'betajuice ~ rank + sexe + age + pjuice ','Upper','interactions','Criterion','bic')

% fitglm(glmtab,'betajuice ~ pjuice ')

glmtab=table;
glmtab.age=zscore(mkidx.age);
glmtab.sexe=zscore(mkidx.sexe);
glmtab.rank=zscore(mkidx.rank);
glmtab.pjuice=zscore(mkidx.pjuice);
glmtab.interc=rg(:,1);
ginter=stepwiseglm(glmtab,'interc ~ rank + sexe + age + pjuice ','Upper','interactions','Criterion','bic')

warning on