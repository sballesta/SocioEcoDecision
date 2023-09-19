%%
function [rez, MDAY,TSTART, MKNAME, MKMOD, MKID, MKNAMEID, CDAY ]= dodays(data,PLOT)
try
    %     clc
    %     keyboard
    tic

    load mkrank_final

    % FIXING NAMES
    data=noname_elo(data,[]);
    mkfinal=noname_elo(mkfinal,[]);

    % spmd
    disp('... transforming data')
    Start=[]; Stop=[]; Duration=[]; Module=[]; Dates=[]; Mname=[];
    % CALCULATING DURATIONS
    for n = 1:size(data,1)
        [h , m, s] = hms(data.StartInstant(n));  Start(n,1) = h*3600 + m*60 + s;
        [h , m, s] = hms(data.EndInstant(n));  Stop(n,1) = h*3600 + m*60 + s;
        Duration(n,1) = Stop(n) - Start(n);
        Module(n,1) = str2num(data.ModuleName{n}(end));
        Dates(n,1) = datenum(data.Date(n));
        Mname{n}=[ data.PlayingSubject{n}];
    end

% if PLOT
%        ConfDate=[]; ConfTime=[]; ConfModule=[]; Cname=[];
%         % CALCULATING DURATIONS
%         for n=1:size(cata,1)
%             h = str2num( cata.Time{n}(1:2) )*3600 ;
%             s = str2num( cata.Time{n}(4:5) )*60 ;
%             ConfTime(n,1)=h+s;
%             ConfDate(n,1) = datenum(cata.Date(n)) ;
%             ConfModule(n,1) =str2num(cata.Module{n}(end)) ;
%             Cname{n}=[ cata.Monkey{n}(1:2) cata.Monkey{n}(end) ' / ' cata.ConflictingMonkey{n}(1:2) cata.ConflictingMonkey{n}(end) ];
%         end
% end

    maxstop=nanmax([Stop ;Start]); maxstart=nanmin(Start);
    msize=numel(maxstart:maxstop);

    disp('... creating days matrix')
    % MDAY=zeros(4,msize);
    % CDAY=zeros(4,msize);
    dates=unique(Dates);
    rez=cell(size(dates,1),2);
    MDAY=cell(size(dates,1),1);
    CDAY=cell(size(dates,1),1);
    TSTART=cell(size(dates,1),1);
    MKNAME=cell(size(dates,1),1);
    MKMOD=cell(size(dates,1),1);
    MKID=cell(size(dates,1),1);
    MKNAMEID=cell(size(dates,1),1);

    for n = 1 : size(dates,1)
        mask = Dates == dates(n);
        tstart = Start(mask); tstop = Stop(mask); tmod = Module(mask); tname=Mname(mask);

        % Creating mkid
        mknameid=zeros(1,numel(tname));
        %         rezname=mkref(1:24,1);
        rezname = mkfinal;
        for mk=1:size(rezname,1)
            %             mask=strcmp(tname,rezname{mk,1});
            mask=contains(tname,rezname{mk,1});
            %         mknameid(mask)=mkrank2(mk,1);
            mknameid(mask)=mk; %mkref{mk,4};
        end
        mknameid(mknameid==0)=NaN;

        % SESS START
        temp=Mname(mask); mksess=cell(4,2);
        for z= 1 : numel(temp)
            mksess{tmod(z),1}=[ mksess{tmod(z),1}; temp{z} ];
            mksess{tmod(z),2}=[ mksess{tmod(z),2}; tstart(z)-maxstart+1 ];
        end

        mday=false(4,msize);
        cday=false(4,msize);
        mkiday=zeros(4,msize);

        for k=1:numel(tstart)
            if ~isnan(tstop(k))
                mday(tmod(k),(tstart(k)-maxstart+1):(tstop(k)-maxstart))=true;
                mkiday(tmod(k),(tstart(k)-maxstart+1):(tstop(k)-maxstart))=mknameid(k);
            end
        end

  
        if PLOT
                            % CONFLICT START
%                 temp=Cname(mask); 
%                 csess=cell(4,2);
%                 for z= 1 : numel(temp)
%                     csess{cmod(z),1}=[ csess{cmod(z),1}; temp{z} ];
%                     csess{cmod(z),2}=[ csess{cmod(z),2}; ctime(z)-maxstart+1 ];
%                 end

            disp([  num2str(n) ' ... plotting: ' datestr(dates(n)) ])
            figure(1)
            %suptitle(datestr(dates(n)))
            
            for k=1:4
                subplot(4,1,k)
                ylabel( ['Module 0' num2str(k)] ); hold on
                plot(mday(k,:),'k');
                plot(cday(k,:).*1.25,'r');

                yidx=1.2:-1/numel(mksess{k,2}):.2;
                yidx=yidx( randperm( numel(mksess{k,2}) ) );
                for z=1:numel(mksess{k,2})
                    %          text(mksess{k,2}(z)-750,rand(1)+.2, mksess{k,1}(z,:)) ;
                    text(mksess{k,2}(z)-750,yidx(z), mksess{k,1}(z,:)) ;
                end

%                 for z=1:numel(csess{k,2})
%                     text(csess{k,2}(z)-50,(rand(1)*.15)+1.15, csess{k,1}(z,:),'Color','r') ;
%                 end

                ylim([0.25 1.5]); xlim([4*3600 22*3600]) ;
            end

            wb = 0;
            while wb~=1
                wb = waitforbuttonpress;
                subtitle(datestr(dates(n)))
            end
            clf
        end
        %%%% FURTHER ANALYSIS
        %%%%% FIND SESSION THAT STARTS ALONE
        %          keyboard
        try
            sessalone=sum(mday);

            if 0 %05/05/2021

                % EXACTLY WHEN START
                sessalone=sessalone(tstart-maxstart+1);
                idx_alone=find(sessalone==1);

            else
                %%% EVENTUALLY ADD MORE STRICT ALONE FILTERING ( e.g ALONE FOR 30 SEC BEFORE STARTING THE SESSION ?)
                %  -n sec before starting
                nsec=20;
                idx=tstart-maxstart-nsec; idx(idx<0)=1;
                %                 idx(idx>=size(sessalone,2))=[]; % FIX
                %                 tstart(idx>=size(sessalone,2))=[];% FIX
                sessalone=sessalone(idx);
                idx_alone=find(sessalone==0);
            end

        catch
            keyboard
        end

        % PACKING DATA
        rez{n,1}=tname(idx_alone)';
        rez{n,2}=tmod(idx_alone);
        %     MDAY=MDAY+mday ;
        %     CDAY=CDAY+cday ;

        MDAY{n}=mday ;
        CDAY{n}=cday ;
        TSTART{n}=tstart-maxstart+1;
        MKNAME{n}=tname;
        MKNAMEID{n}=mknameid;
        MKMOD{n}=tmod';
        MKID{n}=mkiday;

    end

    toc

catch
    keyboard
end
end
