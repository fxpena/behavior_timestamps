classdef BehaviorData
    properties
        box %time stamps
        mouseID %integer or string
        sessionNum %integer
        behavior %label for the appetitive behavior
        identityOfCSplus
        %flavor during the conditioning/test session
        flavor(1,:) char {mustBeMember(flavor,{'chocolate','vanilla'})}='chocolate'; %default is chocolate
        duration    %duration of the behavior session in ms
        graphTitle
        cueDuration=2; %sec
        cueDurationMS=2000; %ms
        trace %sec
        traceMS %ms
        refeeding %volume in mL or quantity in kcal
        volForTest %total volume consumed during consumption test
        weight %grams
        initWeight %initial weight, grams
        feedingFlavor %flavor during the feeding period if there was one
        feedingUnit(1,:) char {mustBeMember(feedingUnit,{'mL','g','kcal'})}='mL'; %default is mL
        %create a struct with the labels and their corresponding number
        dictionary=struct('CS1',1,'CS2',2,'lick',6,'chocolate',7,...
            'vanilla',8,'frame',10,'opto',12);
    end
    properties (Dependent)
        labelOfFlavor
        relativeWeight %that day's weight divided by the initial weight
        feedingCal %amount of feeding in kcal (aka Calories)
    end
    methods
        %constructor function
        function obj = BehaviorData(stamps,mouseID,sessionNum,behavior,CSplus,durationMS,trace)
            %INPUTS
            %stamps: a matrix with 2 columns as defined below
            %    1: (1)Cue1,(2)Cue2,(6)lick,(7)Ensure release,(10)frame
            %    2:timestamp,  ms
            %mouseID: integer or string
            %sessionNum: integer
            %behavior: string
            %CSplus: 'lowTone' or 'highTone'
            %durationMS: duration of the entire behavior session in ms
            %trace: trace interval from CS offset to reward delivery in sec
            if nargin == 0
                % Provide default values so that the object can be constructed
                % without inputs
                obj.box=zeros(2,2);
                obj.mouseID=0;
                obj.sessionNum=0;
                obj.behavior = 'NA';
                obj.identityOfCSplus = 'NA';
                obj.duration=30*60*1000; %duration of full session in ms
            else
                obj.box=stamps;
                if isnumeric(mouseID)
                    obj.mouseID=num2str(mouseID);
                else
                    obj.mouseID=mouseID;
                end
                obj.sessionNum=sessionNum;
                obj.behavior = behavior;
                switch CSplus
                    case {'highTone','high tone','high','tone'}
                        obj.identityOfCSplus ='highTone';
                    case {'lowTone','low tone','low','CS2'}
                        obj.identityOfCSplus = 'lowTone';
                end
                obj.duration=durationMS;
                switch behavior
                    case {'conditioning'}
                        obj.trace=trace; %sec
                    otherwise
                        %for paradigms other than conditioning, automatically
                        %make the trace interval 0
                        trace=0;
                end
                obj.traceMS=trace*1000; %ms
            end
            %Populate additional properties
            obj.graphTitle=['b' obj.mouseID ' ' obj.behavior ' ' num2str(obj.sessionNum)];
            %try to obtain the intial weight from a previous file
            try 
                N=obj.sessionNum-1;
                previous=load(['b' obj.mouseID ' CSPF' num2str(N) '.mat']);
                obj.initWeight=previous.obj.initWeight;
            catch
                %the file couldn't load because it may not exist
            end
            if obj.getNum(8)>0 %if there are vanilla stamps
                obj.flavor='vanilla'; %automatically change the flavor property
            end
        end
        function saveObj(obj)
            %saves the object with a name such as 'b33 CSPF1' and when it
            %loads, it appears as 'obj' in the workspace
            if ischar(obj.mouseID)
                save(['b' obj.mouseID ' CSPF' num2str(obj.sessionNum) '.mat'],'obj')
            else %if it's a number
                ID=num2str(obj.mouseID);
                save(['b' ID ' CSPF' num2str(obj.sessionNum) '.mat'],'obj')
            end
            disp(['data for ' obj.graphTitle ' has been saved'])
        end
        %% Essential methods. Many other methods call these
        function output=getNum(obj,label)
            %Calculate the number of instances for a given event (e.g. #.
            %ensure deliveries)
            %label: (1)Cue1,(2)Cue2,(6)lick,(7)Ensure release,(10)frame
            %   acquisition. Also, 'CS+' and 'CS-' are valid options
            if isstruct(obj.box)
                forCounting=obj.box.regData;
            else
                forCounting=obj.box;
            end
            if isnumeric(label)
                output=sum(forCounting(:,1)==label);
                %if chocolate was asked for, but you meant vanilla
                if label==7 && strcmpi(obj.flavor,'vanilla')
                    output=obj.getNum(8); %vanilla is 8
                end
            else %if the input was a string
                if strcmpi(label,'CS+') || strcmpi(label,'CS-')
                    label=obj.getlabelOfCS(label);
                else %other valid string options
                    label=obj.dictionary.(label);
                end
                output=obj.getNum(label);
            end
        end
        function column=getTimes(obj,label)
            %grab the timestamps of the event of interest (e.g. times of CS
            %onset (in ms)
            %label: (1)Cue1,(2)Cue2,(6)lick,(7)Ensure release,(10)frame
            %acquisition
            if not(isnumeric(label))
                %try to convert the string into it's corresponding number
                label=obj.dictionary.(label);
            end
            indices=obj.box(:,1)==label;
            if any(indices) %if there are stamps with the intended label
                column=obj.box(indices,2);
            else %if there are no stamps with that label
                column=[]; %return an empty matrix
%                 warning(['there are no stamps with the label ' num2str(label) ' for ' obj.graphTitle])
            end
            %if chocolate was asked for, but you meant vanilla
            %even if there was an accidental stamp for chocolate, check for
            %vanilla
            if label==7 && strcmpi(obj.flavor,'vanilla')
                column=obj.getTimes(8); %vanilla is 8
            end
        end
        function label=getlabelOfCS(obj,CS)
            %CS: 'CS+' or 'CS-'
            switch CS
                case 'CS+' %determine the CS+ label
                    switch obj.identityOfCSplus
                        case 'lowTone'
                            label=1;
                        case 'highTone'
                            label=2;
                    end
                case 'CS-' %determine the CS- label
                    switch obj.identityOfCSplus
                        case 'lowTone'
                            label=2;
                        case 'highTone'
                            label=1;
                    end
                otherwise
                    error('whichCS was not properly indicated')
            end
            
        end
        function label=get.labelOfFlavor(obj)
            %Returns the numeric label for the Ensure flavor
            switch obj.flavor
                case 'chocolate'
                    label=7;
                case 'vanilla'
                    label=8;
                otherwise
                    error('obj.flavor is not indicated')
            end
        end
        function w=get.relativeWeight(obj)
            %Returns the ratio of today's weight vs the initial weight
            w=obj.weight/obj.initWeight;
        end
        function k=get.feedingCal(obj)
            %Returns the # of Calories during feeding period
            switch obj.feedingUnit
                %determine caloric density based on units
                case 'mL' %milliliters
                    density=0.93; %kcal/mL
                case 'g' %grams
                    switch obj.feedingFlavor
                        case 'chow'
                            density=3.03; %kcal/g
                        case 'sugar' %sugar pellets
                            density=3.6; %kcal/g
                    end
                case 'kcal' %if it's already in kcal then don't change it
                    density=1; 
            end
            if isempty(obj.refeeding)
                k=0;
            else
                k=density*obj.refeeding;
            end
        end
        function [plus, minus]=getTrialTimes(obj,period,together)
            %Provides the onset and offset times of all trials for both CS's.
            %Input and output in ms.
            %INPUTS
            %period: amount of time after CS onset
            %OPTIONAL INPUT
            %together: set to true to get both CS+ and CS- times in a
            %   single output. They wil be sorted chronologically. the
            %   second output will be empty.
            %OUTPUTS
            %plus: 2 column array for CS+
            %minus: 2 column array for CS-
            
            %check that period is not a vector
            assert(isscalar(period),'input is not a scalar')
            CS1onset=obj.getTimes(1);  %low tone
            CS2onset=obj.getTimes(2);  %high tone
            %if the optional input was given
            if exist('together','var') && together
                %combine the CS times together
                plus=sort([CS1onset; CS2onset]);
                %create 2nd column by adding the period
                plus(:,2)=plus+period; %ms
                %provide the 2nd output
                minus=[];
            else %otherwise proceed with CS+ and CS- separately
                
                %assign CS+ and CS- labels
                switch obj.identityOfCSplus
                    case 'lowTone'
                        CSplusOnset=CS1onset;
                        CSminusOnset=CS2onset;
                    case 'highTone'
                        CSplusOnset=CS2onset;
                        CSminusOnset=CS1onset;
                end
                %calculate the times to close the trials. of course, this
                %closing time is not necessarily the same as the CS offset
                CSplusClose= CSplusOnset+period;
                CSminusClose= CSminusOnset+period;
                %concatenate
                plus=[CSplusOnset CSplusClose];
                minus=[CSminusOnset CSminusClose];
            end
        end
        function [boutOnsets,nLicks]=getBoutTimes(obj,lockout,graph)
            %returns times of the first licks in each of the bouts
            %Approximately 94% or more of licks occur within 2 seconds
            %during conditioning and RT30 so 2 is a good value for the
            %lockout period.
            %INPUTS:
            %lockout: number in sec
            %OPTIONAL INPUT
            %graph: true to show histogram. optional. default is not to graph
            %OUTPUTS
            %boutOnsets: column vector in ms
            %nLicks: # licks in each bout
            
            lickOnsets=obj.getTimes(6); %ms
            totalLickCount=obj.getNum(6);
            ILI=diff(lickOnsets/1000); %inter-lick intervals in sec
            if nargin>2 && graph
                figure; histogram(ILI,25)
                xlabel('inter-lick interval, sec')
                ylabel('count')
                title(obj.graphTitle)
            end
            indices=find(ILI>lockout)+1; %offset by 1
            boutOnsets=lickOnsets(indices); %ms
            boutOnsets=[lickOnsets(1); boutOnsets]; %include the first bout
            %calculate the times of offset -- importantly, this is not equal to the
            %time of the last lick in a bout
            boutOffset=boutOnsets+lockout*1000; %ms
            %group the licks into bouts
            [nLicks,grouped]=BehaviorData.calcEvents(lickOnsets,boutOnsets,boutOffset);
            %if there's a discrepancy between the grouped count and the total
            if sum(nLicks) ~= totalLickCount
                %concatenate the entries of the array into a column
                grouped=vertcat(grouped{:});
                %find out which licks do not fall into the original bins
                indUnaccounted=not(ismember(lickOnsets,grouped));
                %lick onsets that are unaccounted for
                unaccounted=lickOnsets(indUnaccounted);
                moreBtimes=licksIntoBouts(unaccounted,lockout);
                %intregate these times with the other ones
                boutOnsets=sort([boutOnsets;moreBtimes]);
                %recalulate the offset times
                boutOffset=boutOnsets+lockout*1000; %ms
                %recalculate the #licks in each bout
                nLicks=BehaviorData.calcEvents(lickOnsets,boutOnsets,boutOffset);
            end
        end
        function nBouts=getNumBouts(obj,lockout,flag)
            %Counts the number of lick bouts. All licks within a lockout
            %period are considered 1 bout.
            %INPUTS
            %lockout: interval in seconds
            %flag: optional. true to exclude bouts that are too close.
            %   false to include all bouts in the total
            if ~exist('flag','var')
                flag=true;
            end
            bTimes=obj.getBoutTimes(lockout);
            if flag
                bTimes=obj.removeBouts(bTimes,lockout);
            end
            nBouts=length(bTimes);
        end
        function [yes,no]=lickOrNot(obj,whichCS,period,skips)
            %Identifies which trials had a lick during/after CS
            %INPUTS
            %whichCS: 'CS+' or 'CS-'
            %period:time from CS onset, ms
            %skips: optional. a vector of indices for trials that should be 
            %   skipped from this analysis
            %OUTPUTS
            %column vectors of boolean indices
            
            %calculate #licks for all trials
            [nLicks,~,isPlus]=obj.licksAndDel(period); %input in ms
            switch whichCS
                case 'CS+'
                    yes=nLicks(isPlus)>0;
                case 'CS-'
                    yes=nLicks(~isPlus)>0;
                otherwise
                    error('inappropriate input for whichCS')
            end
            
            if nargin >3 && ~isempty(skips) %if skips are provided
                %remove trials that need to be skipped
                yes(skips)=[];
            end
            no=not(yes);
        end
        function v=volPerD(obj)
            %Calculate the average volume for each Ensure delivery
            %OUTPUT
            %v: rate in uL/delivery
            %check that the volume has been added to the object
            assert(~isempty(obj.volForTest),['the volForTest property is empty for ' obj.graphTitle])
            %calculate the volume of the average delivery
            v=obj.volForTest*1000/obj.getNum(obj.labelOfFlavor); %uL/delivery
        end
        %% Analysis and graphing
        function [P,dP,PallTrials,dPallTrials]=rewardProportion(obj,period,override)
            %Calculates the reward proportion experienced during CS+,
            %CS-, and no cue intervals
            %INPUTS
            %period: duration from CS onset to consider for calculation, in
            %   ms. It must be a scalar
            %OPTIONAL INPUT
            %override: By default, if the object contains optogenetic data,
            %   this method subdivides the trials into opto or not. Set
            %   this boolean to true in order to prevent this subdivision.
            %OUTPUTS
            %P: 1x3 vector, mL. column1 is CS+, column2 is CS-, column3 is no cues
            %   if there's optogenetic trials, the order is CS+ w/opto, 
            %   CS+ control, CS- w/opto, CS- control
            %dP: strength of reward contingency. 2 elements, 1 for each CS
            %OPTIONAL OUTPUTS
            %PallTrials: reward probability as a function of trial #. 
            %   struct with 3 fields
            %dPallTrials: reward contingency as a function of trial #.
            %   struct with 2 fields. Contingency is calculated for trial n,
            %   given the reward probability on trials n-1 and n
            % Created by Francisco Pena on 3/5/21
            assert(isscalar(period),'period must be a scalar');
            
            if ~exist('override','var') || isempty(override)
                override=false; %default
            end
            P=zeros(1,3); %row
            %grab the time stamps for Ensure and CS's
            ensureTimes=obj.getTimes(7);
            [plusT,minusT]=obj.getTrialTimes(period);
            %count the # deliveries during each CS presentation
            nRewardsCSplus =BehaviorData.calcEvents(ensureTimes,plusT(:,1),plusT(:,2));
            nRewardsCSminus =BehaviorData.calcEvents(ensureTimes,minusT(:,1),minusT(:,2));
            %collate and sort all CS onset and offset times
            allCS=sortrows([plusT;minusT]);
            %the ITI onsets are the CS offsets
            ITIonset=[0; allCS(:,2)]; %include the time before first CS
            %the ITI offsets are the CS onsets
            ITIoffset=[allCS(:,1); obj.duration]; %include time after last CS
            durations=ITIoffset-ITIonset; %ms
            %identify the intervals that are longer than the CS intervals
            indices=find(durations'>period); %row vector of numeric indices
            %subdivide the ITIs so that they are all the same size
            newIntervals=[];
            for i=indices %for each oversized interval
                %numerate a new row of intervals
                newRow=ITIonset(i):period:ITIoffset(i);
                %pair up the numbers into onsets and offsets
                newColumns=[newRow(1:end-1)' newRow(2:end)'];
                %append the new columns to the larger matrix
                newIntervals=[newIntervals; newColumns];
            end
            %all ITIs that were shorter than period will be omitted
            nRewardsITI=BehaviorData.calcEvents(ensureTimes,newIntervals(:,1),newIntervals(:,2));
            %binarize the trials because we only care about the presence of
            %reward not magnitude
            nRewardsCSplus=nRewardsCSplus>0;
            nRewardsCSminus=nRewardsCSminus>0;
            nRewardsITI=nRewardsITI>0;
            %for a conditioning session
            if strcmp(obj.behavior,'conditioning')
                %modify the counting for the CS+ such that there has to be
                %a lick within 1 sec of reward delivery
                %first calculate the latency to lick after the Ensure deliveries
                latencies = getLatencies(obj,7,6,2000,true);
                nRewardsCSplus=latencies<1;
            end
            %average the booleans, which is equivalent to calculating the
            %proportion of trials with rewards
            P(1)=mean(nRewardsCSplus);
            P(2)=mean(nRewardsCSminus);
            P(3)=mean(nRewardsITI); %baseline
            %calculate the difference relative to baseline
            dP=P(1:2)-P(3);
            %if there was optogenetic stimulation
            if obj.getNum(12)>0 && override==false
                %set the 5th value to baseline for safe-keeping
                P(5)=P(3);
                %subdivide the trials based on whether there was opto
                [yes,no]=obj.optoOrNot('CS+');
                P(1)=mean(nRewardsCSplus(yes)>0); %opto
                P(2)=mean(nRewardsCSplus(no)>0);  %control
                [yes,no]=obj.optoOrNot('CS-');
                P(3)=mean(nRewardsCSminus(yes)>0); %opto
                P(4)=mean(nRewardsCSminus(no)>0);  %control
                %subtract baseline to calculate contingency
                dP=P-P(5);
                dP(5)=[]; %remove this entry
            end
            if nargout>2 %if an additional output was called
                %calculate reward probability as a moving average
                numTrials=1; %# previous trials to use for window
                %only use previous trials and current trial, not forward trials
                %The full window is 2 trials.
                PallTrials.CSplus=movmean(nRewardsCSplus,[numTrials 0]);
                PallTrials.CSminus=movmean(nRewardsCSminus,[numTrials 0]);
                PallTrials.ITI=movmean(nRewardsITI,[numTrials 0]);
                for n=2:length(nRewardsCSplus) %for each CS+ trial
                    %determine the no cue interval that precedes this CS presentation
                    %by using the smallest negative difference from this CS onset time
                    intervalIndex=find(newIntervals(:,1)-plusT(n,1)<0,1,'last'); 
                    %reward probability up to this trial minus ITI
                    dPallTrials.CSplus(n)=PallTrials.CSplus(n)-PallTrials.ITI(intervalIndex);
                end
                %set the first entry to 1
                dPallTrials.CSplus(1)=1;
                for n=2:length(minusT) %for each CS- trial
                    intervalIndex=find(newIntervals(:,1)-minusT(n,1)<0,1,'last');
                    dPallTrials.CSminus(n)=PallTrials.CSminus(n)-PallTrials.ITI(intervalIndex);
                end
                %the first entry is 0 by default and we leave it that way
                %convert both row vectors into columns
                dPallTrials.CSplus=dPallTrials.CSplus';
                dPallTrials.CSminus=dPallTrials.CSminus';
            end
        end
        function [hasEnsure,varargout]=groupBasedOnEnsure(obj,bTimes,lockout)
            %Groups lick bout based on whether there was Ensure present
            %INPUT
            %bTimes: times of lick bout onset, ms
            %lockout: interval that defines the time limit for each bout, sec
            %OUTPUT
            %hasEnsure: boolean. identifies the lick bouts in which there was Ensure
            %       based on the chronological order of the trials
            %OPTIONAL OUTPUT depends on whether type of behavior
            %nDel: a column vector indicating the # deliveries
            %       for each bout
            %delta: for RT30 sessions, this output is the vector of
            %   latencies for lick bouts relative to the matching Ensure 
            %   delivery. in ms
            if isempty(bTimes) %if the bout onsets are not provided
                bTimes=obj.getBoutTimes(lockout); %ms
            end
            bOffset=bTimes+lockout*1000; %ms
            ensureT=obj.getTimes(7); %ms
            nDel=[]; %initialize as empty
            switch obj.behavior
                case 'RT30' %RT30 condition
                %Ensure may have been released before the bout but had not been licked yet
                %or the Ensure may be released during a bout
                    allLickTimes=obj.getTimes(6);
                    delta=NaN(size(ensureT)); %initialize as NaNs
                    keyLicks=[]; %this vector will hold the important lick times
                    for i=1:length(ensureT) %for each delivery
                        %find the first lick after the delivery
                        subtracted=allLickTimes-ensureT(i);
                        if any(subtracted>0) %if there are positive differences
                            %find the smallest one
                            delta(i)=min(subtracted(subtracted>0));
                        end %otherwise, leave this entry as a NaN
                        if ~isnan(delta(i))
                            %add this lick to the list
                            keyLicks=[keyLicks; ensureT(i)+delta(i)];
                        end %o/w leave out this delivery because no lick occurred afterward
                    end
                    %find the bouts that include the important licks
                    hasKeyLick=BehaviorData.calcEvents(keyLicks,bTimes,bOffset);
                    hasEnsure=logical(hasKeyLick);
                    if nargout==2
                    varargout{1}=delta;
                    end
                case {'consumptionTest','consumption test','cued','cuedTest','cued test'}
                    %in FR5 sessions, 2 or more deliveries could occur within 1 bout
                    %include 300ms before bout onset
                    %for all bouts, calculate the #deliveries in each one
                    nDel=BehaviorData.calcEvents(ensureT,bTimes-300,bOffset);
                    %the first delivery is given freely, before any licks, so manually
                    %indicate that the first bout has Ensure
                    nDel(1)=1;
                    hasEnsure=nDel>0;
                    if nargout==2
                    varargout{1}=nDel;
                    end
                otherwise
                    %for non FR5 sessions, the deliveries happen far enough in
                    %time that they cannot occur within 1 lick bout. And we
                    %continue the assumption that if a drop is not caught
                    %immediately, the mouse will not receive it
                    
                    %for all bouts, calculate the #deliveries in each one
                    nDel=BehaviorData.calcEvents(ensureT,bTimes,bOffset);
                    hasEnsure=nDel>0;
                    if nargout==2
                    varargout{1}=nDel;
                    end
            end
        end
        function E=getExtraLicks(obj)
            %Calculate the # of non-eating licks
            %OUTPUT:
            %E: the # of extra licks
            %NOTE: this method may not be appropriate for non RT sessions
            totalNumLicks=obj.getNum(6);
            lockout=2;
            %calculate the # licks in bouts
            [bTimes,boutLicks]=obj.getBoutTimes(lockout); %2 second lockout
            %boutLicks is a column vector
            %determine which bouts have Ensure and/or are the first bouts after a delivery
            hasEnsure=groupBasedOnEnsure(obj,bTimes,lockout);
            %subtract the lick counts of the bouts that have Ensure
            E=totalNumLicks-sum(boutLicks(hasEnsure)); 
        end
        function [newBtimes,ind]=removeBouts(obj,bTimes,lockout,thresh,graph)
            %Remove bout times that occurr within 0.3 sec of a previous bout
            %INPUT
            %bTimes: should be all the bouts. a vector of times, ms
            %lockout: in sec
            %thresh: required amount of time b/w bouts, in ms
            %graph: optional. true to plot a histogram of interbout times
            %OUTPUT
            %newBtimes:vector of bout times that already excludes the
            %   undesired bouts. ms
            %ind: numeric indices that can be plugged into the original
            %   list of bouts to produce newBtimes
            %NOTE: Previously, this method removed times that occurr within 
            %0.3 sec of a previous bout. on 7/14/21 I changed it to remove
            %based on the lockout interval
            if ~exist('thresh','var') || isempty(thresh)
                %default is to use the lockout period
                thresh=lockout*1000; %convert s to ms
            end
            if ~exist('graph','var')
                graph=false; %default
            end
            if isempty(bTimes) %if the bout times are not provided
                bTimes=obj.getBoutTimes(lockout);
            end
            %what's the duration b/w each offset and the next bout?
            interBout=diff(bTimes)-lockout*1000; %ms
            %only keep the bouts that have sufficient separation
            b=[true; interBout>thresh]; %also include the first bout by default
            newBtimes=bTimes(b); %index into the original times
            if nargout==2
                ind=find(b); %convert to numeric indices
            end
            if graph
            figure; histogram(interBout/1000,'BinEdges',[0:0.5:10 inf])
            xlabel('inter-bout interval, sec')
            ylabel('# bouts')
            title(obj.graphTitle)
            end
        end     
        function hasCS=groupBasedOnCS(obj,bTimes,interval,before)
            %Group lick bouts based on whether there was a CS+
            %before or during
            %INPUT all in ms
            %bTimes: times of lick bout onset, ms
            %interval: lockout period using to calculate lick bouts, ms
            %OUTPUT
            %hasCS: indices of lick bouts in which there was CS+
            %       based on the chronological order of the bouts
            label=obj.getlabelOfCS('CS+');
            plusT=obj.getTimes(label);
            limitFirst=bTimes-before; %first time boundary for all bouts
            limitSecond=bTimes+interval; %second boundary for all bouts
            %determine how many CS onsets are within the time limits for all bouts
            numCS =BehaviorData.calcEvents(plusT,limitFirst,limitSecond);
            hasCS=find(numCS>0);
        end
        function showBoutCount(obj)
            %Plots the #lick bouts that have a CS+ beforehand
            lockout=2; %sec
            [bTimes,~]=obj.getBoutTimes(lockout); %ms
            durations=(1:10)*1000; %ms
            nBouts=zeros(size(durations)); %initialize
            for i=1:length(durations)
                indices=obj.groupBasedOnCS(bTimes,0,durations(i));
                %calculate the #bouts based on this duration
                nBouts(i)=length(indices);
            end
            durations=durations/1000; %convert to sec
            figure; scatter(durations,nBouts,'filled')
            xlabel('pre-interval, sec'); xlim([0 max(durations)])
            ylabel('# lick bouts'); ylim([0 inf])
            title(obj.graphTitle)
        end
        function deltas=nearestLick(obj,cutoff,graph)
            %Calculates the durations between deliveries and first licks
            %OPTIONAL INPUT
            %graph: true to plot histogram
            %cutoff: duration in sec. latencies longer than this value will
            %   be removed from the output
            %OUTPUT
            %deltas: vector in sec that lists the latency to lick for each
            %   delivery of Ensure
            if ~exist('graph','var') || isempty(graph)
                graph=false; %default
            end
            if ~exist('cutoff','var') || isempty(cutoff)
                cutoff=180; %the default cut off is 180 sec
            end
            
            allLickTimes=obj.getTimes(6);
            ensureT=obj.getTimes(7);
            deltas=[]; %this vector will hold the differences
            for i=1:length(ensureT) %for each delivery
                %find the first lick after the delivery
                subtracted=allLickTimes-ensureT(i); %make this delivery the reference
                delta=min(subtracted(subtracted>0)); %smallest positive difference
                if ~isempty(delta)
                    deltas=[deltas; delta]; %add this diff to the list
                end %o/w leave out this delivery because no lick occurred afterward
            end
            deltas=deltas/1000; %convert from ms to sec
            %omit outliers and/or cases where a lick was not counted
            deltas(deltas>cutoff)=[];
            if graph
                %1 second bins up to 15 sec. 1 bin for all other durations
                binEdges=[0:15 inf];
                figure;
                histogram(deltas,binEdges)
                xlabel('time b/w delivery and lick, sec')
                xlim([0 15]); ylim([0 20])
                ylabel('# trials')
                title(obj.graphTitle)
            end
        end
        function latencies = getLatencies(obj,label,event,cap,overrideRem)
            %Calculates the latency of an event after the given cue
            %presentations in sec
            %INPUT
            %label: 1 or 2 (6 and 7 are also possible). or a vector of
            %   timestamps in ms
            %event: event of interest. 6 for licks, 7 for Ensure release
            %cap: the maximum duration to allow for latencies such as the cue
            %       duration, in ms
            %overrideRem: optional. by default, latencies longer than cap 
            %   are removed from the output vector. Set to true to override
            %   the removal of such latencies
            %OUTPUT
            %latencies:column vector of times, in sec
            if ~exist('overrideRem','var') %if it wasn't provided
                overrideRem=false; %default is to remove the latencies that go over the cap
            end
            %if the obj is a dummy
            if strcmp(obj.behavior,'NA')
                latencies=[]; %return an empty output
            else %otherwise proceed with the function
                if length(label)==1
                    %collect the time stamps of interest
                    referenceTimes=obj.getTimes(label);
                else %if it's a vector
                    referenceTimes=label; %pass it along
                end
                eventTimes=obj.getTimes(event);
                
                if ~isempty(referenceTimes) %if there are onsets for the desired event
                    numCue=length(referenceTimes);
                    latencies=zeros(numCue,1); %initialize
                    for i=1:numCue %for each cue
                        allDiff=eventTimes-referenceTimes(i); %times relative to 1 cue
                        try
                            %rule out the negative differences
                            latencies(i)=min(allDiff(allDiff>0));
                        catch
%                             warning('There may not have been a lick for a cue presentation')
                            latencies(i)=cap;
                        end
                        
                    end
                    latencies(latencies>cap)=cap; %cap the latencies that went over
                    if not(overrideRem)
%                     if event==7 || event==8 %for Ensure events
                        %remove the "latencies" for events that didn't occur or went over
                        %because those events wouldn't make sense
                        latencies(latencies==cap)=[];
%                     end
                    end
                else
                    latencies=[];
                    disp(['Cue/event ' num2str(label) ' was not presented during this session.'])
                end
                %convert ms to sec
                latencies=latencies/1000;
            end
        end
        function h=plotLatencies(obj,event,latencies)
            %Plots latencies and saves the 2 panel graph
            %INPUT
            %event: 6, 7, or 8
            %latencies: in sec. output from getLatencies
            %OUTPUT
            %h: handle for the figure
            h=figure; subplot(2,1,1)
            histogram(latencies,10) %10 bins
            switch event
                case 6
                    xlabel('lick latency, sec')
                    ending='lick latencies';
                case {7,8} %chocolate or vanilla
                    xlabel('ensure latency, sec')
                    ending='ensure latencies';
            end
            ylabel('# trials')
            title(obj.graphTitle)
            subplot(2,1,2)
            ecdf(latencies)
%             savefig(h, [obj.graphTitle ' ' ending]);
        end
        function output=meanLatency(obj,numTrials)
            %Calculate the average latency to lick after the CS+. Only uses
            %the last X trials (e.g. last 10 trials)
            % numTrials: number of trials to use for the calculation
            %             write 0 to include all trials
            % output: in sec
            
            %consider moving this switch-case block to the constructor function
            switch obj.identityOfCSplus
                case 'lowTone'
                    label=1;
                case 'highTone'
                    label=2;
            end
            %cap the latencies at 120 sec because that is the maximum ITI
            latencies = obj.getLatencies(label,6,120000);
            totalNumTrials=obj.getNum(label);
            if numTrials==0
                %use all trials
                output=mean(latencies);
            elseif numTrials >0 && numTrials < totalNumTrials
                %if it is the first session
                if obj.sessionNum==1
                    %use the first x trials
                    output=mean(latencies(1:numTrials));
                else %otherwise
                    %use the last x trails of the session
                    %add 1 to the index to avoid including an additional trial
                    index=totalNumTrials-numTrials+1;
                    output=mean(latencies(index:end));
                end
            else
                error('numTrials is out of bounds')
            end
            
        end
        
        %% Graphing raw data (mostly)
        function plotWholeSession(obj,useNew,showBouts)
            %Show all CS presentaions, licks, and Ensure deliveries
            %OPTIONAL INPUTS
            %useNew: true to use a new figure. false to use the current fig
            %omitBouts: true to include asterisks for lick bouts
            if ~exist('useNew','var') || isempty(useNew)
                useNew=true; %default
            end
            if ~exist('omitBouts','var') || isempty(showBouts)
                showBouts=false; %default
            end
            %calculate the onsets and offsets of all CS presentations
            [plus, minus]=obj.getTrialTimes(obj.cueDurationMS); %ouput is in ms
            plus=plus/1000; minus=minus/1000; %convert ms to sec
            yMax=1;
            Y=[-1/15 -1/15]*yMax; %y values for all CS lines
            if useNew
            figure; hold on
            xlabel('sec')
            ylim([-0.1 inf])
            end
            if ~isempty(plus) %if there are CS+
            %plot a line for the duration of each CS+ presentation
            plot(plus,Y,'color',[0 0.7 0.6],'Linewidth',3) %teal
            end
            if ~isempty(minus) %if there are CS-
            %plot a line for the duration of each CS- presentation
            plot(minus,Y,'color',[0.1 0.1 0.1],'Linewidth',3) %yellow
            end
%             %licks in black ticks
%             plot(obj.getTimes(6)/1000,zeros(obj.getNum(6),1),'k+')
            if showBouts
                %calculate the lick bout onset times
                lockout=2; %2sec lockout
                bTimes=obj.getBoutTimes(lockout); 
                nBouts=length(bTimes);
                has=obj.groupBasedOnEnsure(bTimes,lockout); %boolean indices
                %lick bout onsets
                plot(bTimes/1000,ones(nBouts,1),'k+')
                %blue for the ones that have Ensure
                plot(bTimes(has)/1000,ones(sum(has),1),'b*') 
            end
            %Ensure delivery in magenta ticks
            plot(obj.getTimes(7)/1000,zeros(obj.getNum(7),1),'m+')
            %calculate the reward rate with a sliding window
            window=8; %interval in seconds
            [rollin,t]=obj.rollingCount(7,window);
            plot(t,rollin)
            ylabel('# rewards per 8 sec')
            xlim([0 obj.duration/1000])
            title(obj.graphTitle)
        end
        function plotWholeSessionGUI(obj,ax)
            %Show all CS presentaions, licks, and Ensure deliveries
            %INPUT
            %ax: axes object
            
            %calculate the onsets and offsets of all CS presentations
            [plus, minus]=obj.getTrialTimes(obj.cueDurationMS); %ouput is in ms
            plus=plus/1000; minus=minus/1000; %convert ms to sec
            yMax=1;
            Y=[-1/15 -1/15]*yMax; %y values for all CS lines
            if ~isempty(plus) %if there are CS+
            %plot a line for the duration of each CS+ presentation
            plot(ax,plus,Y,'color',[0 0.7 0.6],'Linewidth',3) %teal
            end
            if ~isempty(minus) %if there are CS-
            %plot a line for the duration of each CS- presentation
            plot(ax,minus,Y,'color',[0.9 0.9 0],'Linewidth',3) %yellow
            end
            %calculate the lick bout onset times
            lockout=2; %2sec lockout
            bTimes=obj.getBoutTimes(lockout); 
            nBouts=length(bTimes);
            h=obj.groupBasedOnEnsure(bTimes,lockout);
            %licks in black ticks
            plot(ax,obj.getTimes(6)/1000,zeros(obj.getNum(6),1),'k+')
            %lick bout onsets
            plot(ax,bTimes/1000,repmat(Y(1),nBouts,1),'b+','LineWidth',2)
            plot(ax,bTimes(h)/1000,repmat(Y(1),sum(h),1),'b*','LineWidth',2)
            %Ensure delivery in magenta ticks
            plot(ax,obj.getTimes(7)/1000,zeros(obj.getNum(7),1),'m+')
            
        end
        function plotLickPSTHsmooth(obj,binWidth,event,yMax)
            %Plots the smoothed lick PSTH of CS+ and CS- or just aligned to
            %lick bouts
            %INPUTS
            %binWidth: width in sec for sliding window
            %event: 'both CS' or 'lick bout'
            %yMax: value for setting the height of the vertical lines at CS
            %   onset and offset
            switch event
                case {'both CS','both', 'CS'}
                    %calculate the PSTH for each CS
                    [psthPlus,t,allTrialsP]=obj.lickPSTHsmooth('CS+',binWidth);
                    [psthMinus,~,allTrialsM]=obj.lickPSTHsmooth('CS-',binWidth);
                    %plot against time
                    figure; hold on
                    %CS+ in teal
                    plot(t,psthPlus,'Color','k','LineWidth',2);
                    %CS- in yellow
                    plot(t,psthMinus,'Color','k','LineWidth',2);
                    shadeSEM(t',allTrialsP,[0 0.7 0.6]);
                    shadeSEM(t',allTrialsM,[0.9 0.9 0]);
                    %dashed lines at CS onset and offset
                    plot([0 0],[0 yMax],'k--')
                    plot([2 2],[0 yMax],'k--')
                    legend('CS+','CS-','Location','Best')
                    xlabel('time, s')
                    ylabel('licks/s')
                    title(obj.graphTitle)
                case 'CS+'
                    %calculate the PSTH
                    [psthPlus,t,allTrials]=obj.lickPSTHsmooth('CS+',binWidth);
                    %plot against time
                    figure; hold on
                    %CS+ in teal
                    plot(t,psthPlus,'Color',[0 0.7 0.6],'LineWidth',2);
                    shadeSEM(t',allTrials,[0 0.7 0.6]);
                    %dashed lines at CS onset and offset
                    plot([0 0],[0 yMax],'k--')
                    plot([2 2],[0 yMax],'k--')
                    xlabel('time, s')
                    ylabel('licks/s')
                    title(obj.graphTitle)
                case 'CS-'
                    %calculate the PSTH
                    [psth,t,allTrials]=obj.lickPSTHsmooth('CS-',binWidth);
                    %plot against time
                    figure; hold on
                    %CS+ in teal
                    plot(t,psth,'Color','k','LineWidth',2);
                    shadeSEM(t',allTrials,[0.9 0.9 0]);
                    %dashed lines at CS onset and offset
                    plot([0 0],[0 yMax],'k--')
                    plot([2 2],[0 yMax],'k--')
                    xlabel('time, s')
                    ylabel('licks/s')
                    title(obj.graphTitle)
                case {'lick bout','bout'}
                    %calculate the lick PSTH
                    [lickRate,t]=obj.lickPSTHsmooth('lick bout',binWidth,[],[]);
                    %plot against time
                    figure; hold on
                    plot(t,lickRate,'k','LineWidth',2)
                    %dashed line at first lick
                    plot([0 0],[0 5],'k--')
                    xlabel('time, s')
                    ylabel('licks/s')
                    title(obj.graphTitle)
                otherwise
                    error('Inappropriate input for event')
            end
        end
        function [PSTH,t,lickRates]=lickPSTHsmooth(obj,event,binWidth,before,after)
            %Compute a lick PSTH with a sliding window
            %INPUTS
            %event: 'CS+', 'CS-', or 'lick bout'
            %binWidth: width of the moving window in sec
            %before, after: # seconds. is empty vectors to use default
            %   values from genRaster(). If you call this method without
            %   these inputs at all, the defaults are set by this method
            %   instead of genRaster()
            %OUTPUTS
            %PSTH: row vector. the unit is licks/s
            %t: vector of time points for each value in PSTH. s, not ms
            %OPTIONAL OUTPUT
            %lickRates: matrix. each row is the moving average of a trial
            if ~exist('before','var')
                before=6; %default
            end
            if ~exist('after','var')
                after=12; %default
            end
            binsPerSec=10; %10 bins/s so that the step size is 100ms
            %indicate the # seconds before and after CS onset or leave
            %those inputs empty to use the default [-2 10]
            [~,lickTicks,~,t]=obj.genRaster(event,before,after);
            %convert t from ms to s
            t=t/1000;
            %add the durations before and after CS onset to get the total
            numSec=abs(t(1))+t(end);
            %lickTicks has a time bin at each ms. group the ticks into 
            %non-overlapping bins of 100ms
            lickTicks=binMat(lickTicks,binsPerSec*numSec); % 10 bins/s
            %downsample the time vector as well, 100 ms increments
            t=t(1:100:end-1); %omit the last bin b/c binMat often cuts off the last time bin
            %first average across trials
            PSTH=mean(lickTicks); %average across rows (along columns)
            %take the moving average so that each value is the average lick
            %rate for the set of bins centered at that bin
            nBins=binsPerSec*binWidth;
            PSTH=movmean(PSTH,nBins)*1000; 
            %the values are multiplied by 1000 at the end to bring the
            %units to licks/s. I have used the calcLickRates() method to
            %verify that these rates make sense empirically.
            if nargout==3
                %take the moving average of each trial
                lickRates=movmean(lickTicks,nBins,2)*1000; %move along the rows
            end
        end
        function plotLickPSTH(obj,binWidth,bothCS)
            %Plots the lick PSTH for CS+ and CS-
            %INPUTS
            %binWidth: width for each bin in sec
            %bothCS: true to include CS- with CS+
            if ~exist('bothCS','var')
                bothCS=false; %default
            end
            [psthPlus,edges]=obj.lickPSTH('CS+',binWidth);
            [psthMinus,~]=obj.lickPSTH('CS-',binWidth);
            edges=edges/1000; %convert from ms to sec
            h=figure; hold on
            histogram('BinEdges',edges,'BinCounts',psthPlus,'FaceColor',[0 0.7 0.6]); %CS+ in teal
            if bothCS
            histogram('BinEdges',edges,'BinCounts',psthMinus,'FaceColor',[0.9 0.9 0]); %CS- in yellow
            legend('CS+','CS-','Location','Best')
            end
            ylabel('# licks/trial'); ylim([0 4])
            xlabel('time, sec'); xlim([edges(1) edges(end)])
            
            currentTitle=[obj.graphTitle ' lick PSTH'];
            title(currentTitle)
        end
        function [PSTH,edges]=lickPSTH(obj,whichCS,binWidth)
            %Calculates the lick peri-stimulus time histogram
            %INPUTS
            %whichCS: 'CS+' or 'CS-'
            %binWidth: width for each bin in sec
            %OUTPUTS
            %PSTH: row vector
            %edges: for the bins of the histogram
            binWidth=binWidth*1000; %convert from sec to ms
            [~,lickTicks,~,time]=obj.genRaster(whichCS);
            %only use the trials where there were licks during/after CS
            [yes,~]=obj.lickOrNot(whichCS,8000); %8000 ms
            lickTicks=lickTicks(yes,:);
            nTrials=size(lickTicks,1);
            %time is a vector of timestamps relative to CS onset
            %Turn it into a matrix by repeating it
            T=repmat(time,nTrials,1); %1 row for each trial
            %use the 1s in the tick matrix to create a timestamp matrix
            %in other words, filter out all the timestamps that do not have
            %a lick
            lickTimes=lickTicks.*T; %element-wise product
            %collapse into a single column and remove the zeros
            lickTimes=lickTimes(:); lickTimes(lickTimes==0)=[];
            %bin the data based on the window
            edges=time(1):binWidth:time(end); %edges of the bins
            [N,edges]=histcounts(lickTimes,edges); %sum across trials
            PSTH=N/nTrials; %average across trials
        end
        function [withB,withoutB,t]=lickPSTHbout(obj,whichCS,binWidth)
            %Compute a lick PSTH with a sliding window
            %INPUTS
            %whichCS: 'CS+' or 'CS-'
            %binWidth: width of the moving window in sec
            %OUTPUTS
            %withB: lick rate over time, licks/s. only for trials with
            %   a lick bout within 8 s of CS onset
            %withoutB: only for trials without bouts after
            %t: vector of time points for each value in PSTH. s, not ms
            
            binsPerSec=10; %10 bins/s so that the step size is 100ms
            %indicate the # seconds before and after CS onset or leave
            %those inputs empty to use the default [-2 10]
            [~,lickTicks,~,t]=obj.genRaster(whichCS,6,12);
            %convert t from ms to s
            t=t/1000;
            %add the durations before and after CS onset to get the total
            numSec=abs(t(1))+t(end);
            %lickTicks has a time bin at each ms. group the ticks into 
            %non-overlapping bins of 100ms
            lickTicks=binMat(lickTicks,binsPerSec*numSec); % 10 bins/s
            %downsample the time vector as well, 100 ms increments
            t=t(1:100:end-1); %omit the last bin b/c binMat often cuts off the last time bin
            %split up the trials based on licks
            [yes,no]=lickOrNot(obj,whichCS,8000);
            %average across trials using the boolean indices
            withB=mean(lickTicks(yes,:)); %average across rows (along columns)
            withoutB=mean(lickTicks(no,:));
            %take the moving average so that each value is the average lick
            %rate for the set of bins centered at that bin
            nBins=binsPerSec*binWidth;
            withB=movmean(withB,nBins)*1000;
            withoutB=movmean(withoutB,nBins)*1000;
            %the values are multiplied by 1000 at the end to bring the
            %units to licks/s. I have used the calcLickRates() method to
            %verify that these rates make sense empirically.
        end
        function p=compareBaselineLicks(obj,interval,option)
            %Compare the # licks in a given interval relative to CS onset
            %INPUT
            %whichCS: 'CS+' or 'CS-'
            %interval: 2-element vector e.g. [-4 -1] in sec
            %option: this string determines which types of trials to compare
            %   'licks' compares those with licks vs those without
            %   'CS' compares the CS+ trials vs CS- trials
            %OUTPUT
            %p: p value of a t-test. Might be a struct with 1 value for CS+
            %   and 1 value for CS-
            %NOTE: comparing the # of licks is equivalent to comparing the
            %lick rates because both intervals have the same duration
            
            intForLick=8000; %ms
            lickTimes=obj.getTimes(6); %grab the lick timestamps
            interval=interval*1000; %convert from s to ms
            %grab the CS onsets and add the 2nd element of the interval
            [plusT,minusT]=obj.getTrialTimes(interval(2));
            %alter the onset times by the 1st element of the interval
            plusT(:,1)=plusT(:,1)+interval(1);
            minusT(:,1)=minusT(:,1)+interval(1);
            %count the # licks in each of these intervals
            N.plus=BehaviorData.calcEvents(lickTimes,plusT(:,1),plusT(:,2));
            N.minus=BehaviorData.calcEvents(lickTimes,minusT(:,1),minusT(:,2));
            switch option
                case {'licks','lick bouts','bouts'}
                    %seperate the trials based on having licks or not
                    [yes,no]=obj.lickOrNot('CS+',intForLick);
                    %compute the probability
                    [~,p.plus]=ttest2(N.plus(yes),N.plus(no));
                    [yes,no]=obj.lickOrNot('CS-',intForLick);
                    [~,p.minus]=ttest2(N.minus(yes),N.minus(no));
                case 'CS'
                    %compare the CS+ and CS- trials
                    [~,p]=ttest2(N.plus,N.minus);
                otherwise
                    error('Inappropriate input for option')
            end
        end      
        function [lickTicks,ensureTicks]=createTicks(obj,start,stop,nTrials)
            %Creates matrices that are the ticks in a raster plot
            %INPUTS
            %start, stop: column vectors that define how to carve the
            %   trials, ms
            %nTrials: optional. it can be lower than the full number of trials
            %NOTE: this function has essentially the same algorithm as
            %calcEvents()
            if ~exist('nTrials','var') || isempty(nTrials)
                nTrials=length(start); %default is to use all trials
            end
            %grab the times of the events of interest
            lickTimes=obj.getTimes(6);
            ensureTimes=obj.getTimes(7);
            %each row is a trial across time
            lickTicks=zeros(nTrials,length(start(1):stop(1)));
            ensureTicks=lickTicks; %same initialization as above
            for i=1:nTrials %for each trial
                %these times must be within the time limits of the whole session
                if start(i)>0 && stop(i)<obj.duration
                    %create a dummy vector with a timepoint for every ms
                    allTimesInWindow=start(i):stop(i);
                    %place a 1 at every time with a lick
                    lickTicks(i,:)=ismember(allTimesInWindow,lickTimes);
                    %place a 1 at every time with ensure release
                    ensureTicks(i,:)=ismember(allTimesInWindow,ensureTimes);
                end
            end
        end
        function [nTrials,lickTicks,ensureTicks,varargout]=genRaster(obj,event,before,after,trialLim)
            %Generates the data for raster plots of licks and ensure
            %deliveries, but the graph is created by another function
            %INPUTS
            %event: 'CS+','CS-', 'lick bout', or 'Ensure'
            %OPTIONAL INPUT
            %before: # of seconds before CS onset or event onset
            %after: # of seconds after CS onset
            %if before and after are not provided, default values are
            %determined based on the event
            %trialLim: if this isn't provided, all trials are used
            %          if this is an integer, n, the first n trials are used
            %          if this is a 2 element matrix,[n1,n2], trials n1 to n2
            %          are used for the plot
            %OUTPUTS
            %nTrials: #trials included
            %lickTicks: matrix where each row is a trial over time. mostly
            %   zeros, but there is a 1 at every ms where there is a lick
            %ensureTicks: same as lickTicks but for Ensure deliveries
            %time:row vector from -2,000 to 10,000 ms. optional output
            switch event
                case {'CS+','CS-'}
                    %grab event times
                    CSonset=obj.getTimes(obj.getlabelOfCS(event));
                    
                    %if trialLim was provided
                    if exist('trialLim','var') && ~isempty(trialLim) 
                        %set upper limit based on the input
                        nTrials=min(length(CSonset),max(trialLim));
                    else
                        %count # of trials (e.g. # CS presentations)
                        nTrials=length(CSonset); %use all trials
                    end
                    %if before was not provided
                    if ~exist('before','var') || isempty(before)
                        before=2000; %default in ms
                    else
                        %convert input from s to ms
                        before=before*1000;
                    end
                    %if after was not provided
                    if ~exist('after','var') || isempty(after)
                        after=10000; %default in ms
                    else
                        after=after*1000; %ms
                    end
                    time= -(before):after;
                    %determine the start and stop times of the epochs
                    start=CSonset-before;
                    stop=CSonset+after;
                case 'lick bout'
                    before=2000; after=6000; %ms
                    lockout=2; %2 sec
                    bOnsets=obj.getBoutTimes(lockout);
                    %filter out the bouts that occur too soon after another
                    bOnsets=obj.removeBouts(bOnsets,lockout);
                    start=bOnsets-before;
                    stop=bOnsets+after;
                    nTrials=length(start);
                    time= -(before):after;
                case 'lick bout Ensure'
                    before=4000; after=6000; %ms
                    lockout=2; %2 sec
                    bOnsets=obj.getBoutTimes(lockout);
                    %filter out the bouts that occur too soon after another
                    [~,allowedInd]=obj.removeBouts(bOnsets,lockout);
                    indGroup=obj.groupBasedOnEnsure(bOnsets,lockout);
                    indGroup=find(indGroup); %convert from boolean to numeric
                    %use allowedInd to create a conjuction of both the grouping
                    %and being allowed
                    both=indGroup(ismember(indGroup,allowedInd));
                    withEns=bOnsets(both);
                    start=withEns-before;
                    stop=withEns+after;
                    nTrials=length(start);
                    time= -(before):after;
                case {'Ensure','ensure','ens'}
                    before=2000; after=6000; %ms
                    %grab all the delivery times
                    T=obj.getTimes(7);
                    start=T-before;
                    stop=T+after;
                    nTrials=length(start);
                    time= -(before):after;
                otherwise
                    error('inappropriate input for event')
            end
            [lickTicks,ensureTicks]=obj.createTicks(start,stop,nTrials);
            if nargout==4
                varargout{1}=time;
            end
        end
        function rasterPlot(obj,event,sorting,trialLim)
            %Creates a raster plot of licks aligned to the desired event. 
            %Everything in the function is in ms but the graph will have 
            %units of sec.
            %INPUTS
            %event: 'CS+','CS-', 'lick bout', or 'Ensure'
            %OPTIONAL INPUTS
            %trialLim: if this isn't provided, all trials are used
            %           if this is an integer, n, the first n trials are used
            %           if this is a 2 element matrix,[n1,n2], trials n1 to n2
            %           are used for the plot
            %sorting: allows you to sort the CS trials based on whether 
            %   there are licks during/after the CS or whether there was 
            %   optogenetic stimulation
            % Created by Francisco Pena, 8/22/18
            if ~exist('sorting','var')
                sorting='default'; %default is to plot trials chronologically
            end
            
            %do not proceed with function if the obj is a dummy
            assert(~strcmp(obj.behavior,'NA'),'BehaviorData object is a dummy');
            if exist('trialLim','var') && ~isempty(trialLim) %if trialLim was provided
                [nTrials,lickTicks,ensureTicks,time]=obj.genRaster(event,trialLim);
            else
                [nTrials,lickTicks,ensureTicks,time]=obj.genRaster(event);
            end
            
            h=figure;  hold on
            %sort the CS trials if desired
            switch sorting
                %NOTE: sorting only makes sense for if event is CS+ or CS-
                case {'opto','optogenetics'}
                    [yes,no]=optoOrNot(obj,event);
                    %stack the indices
                    orderedInd=[find(no); find(yes)]; %numeric indices
                    %change order of ticks
                    lickTicks=lickTicks(orderedInd,:);
                    ensureTicks=ensureTicks(orderedInd,:);
                    %plot a dashed line to divide these trials
                    plot(time, repmat(sum(no),size(time)),'k--')
                case {'licks','lick bouts','bouts','lick bout'}
                    [yes,no]=lickOrNot(obj,event,8000); %8000 ms
                    %stack the indices
                    orderedInd=[find(no); find(yes)]; %numeric indices
                    %change order of ticks
                    lickTicks=lickTicks(orderedInd,:);
                    ensureTicks=ensureTicks(orderedInd,:);
                    %plot a dashed line to divide these trials
                    plot(time, repmat(sum(no),size(time)),'k--')
            end
            %plot all trials stacked on each other
            switch event
                case {'CS+','CS-'}
                    %draw a rectangle for the duration of the CS
                    xRect=[0 obj.cueDurationMS obj.cueDurationMS 0]/1000; %sec
                    yRect=[0 0 nTrials+0.5 nTrials+0.5];
                    switch event
                        case 'CS+'
                            %create a handle for the rectangle
                            hRect=fill(xRect,yRect,[0 0.7 0.6]); %CS+ in teal
                            currentTitle=[obj.graphTitle ' CS+'];
                        case 'CS-'
                            hRect=fill(xRect,yRect,[0.9 0.9 0]); %CS- in yellow
                            currentTitle=[obj.graphTitle ' CS-'];
                    end
                    hRect.FaceAlpha=0.5; %set the transparency to 50%
                    hRect.EdgeColor='none'; %remove the outline
                case 'lick bout'
                    currentTitle=[obj.graphTitle ' lick bouts'];
                case {'Ensure','ensure','ens'}
                    currentTitle=[obj.graphTitle ' Ensure'];
                otherwise
                    currentTitle=obj.graphTitle;
            end
            %the X stays the same and the Y is the trial number
            Y=diag(1:nTrials)*lickTicks;
            Yensure=diag(1:nTrials)*ensureTicks;
            %make all the zero entries negative so they don't show on the plot
            Y(Y==0)=-2; Yensure(Yensure==0)=-2;
            time=time/1000; %convert from ms to sec
            plot(time,Y,'+k') %licks in black
            plot(time,Yensure,'*m') %ensure in magenta
            ylabel('trial');
            %change the amount of trials visibly plotted based on input
            if nargin==5 && ~isempty(trialLim) %if trialLim was provided
                %check whether trialLim has 2 numbers
                if length(trialLim)==2
                    %sort is here in case the numbers are out of order
                    ylim(sort(trialLim))
                else %if there's only 1 value provided
                    ylim([1 trialLim])
                end
            else %otherwise show all the trials
                ylim([1 nTrials])
            end
            xlabel('time, sec'); 
            xlim([time(1) time(end)])
            title(currentTitle)
        end        
        %% Analysis of cued consumption tests
        function [v,varargout]=volDuring(obj,start,stop)
            %Calculate the volume consumed during a specific window of time
            %INPUTS
            %start, stop: times in sec.
            %OUTPUTS
            %v: volume in mL
            %nDel: optional. # deliveries during this window
            
            volumePerDelivery=volPerD(obj)/1000; %convert from uL to mL
            ensureTimes=obj.getTimes(7)/1000; %convert ms to sec
            %count how many deliveries occurred during the window
            nDel=BehaviorData.calcEvents(ensureTimes,start,stop); %all inputs have the same time units
            %multiply #deliveries by delivery rate
            v=nDel*volumePerDelivery;
            if nargout==2
                varargout{1}=nDel;
            end
        end
        function [rate,varargout]=lickDuring(obj,start,stop)
            %Calculate the lick rate during a specific window of time
            %INPUTS
            %start, stop: times in sec.
            %OUTPUTS
            %rate: lick rate during the desired time, lick/s
            %nLicks: optional. # licks during that time
            
            lickTimes=obj.getTimes(6)/1000; %convert ms to sec
            %count how many licks occurred during the window
            nLicks=BehaviorData.calcEvents(lickTimes,start,stop); %all inputs have the same time units
            rate=nLicks/(stop-start);
            if nargout==2
                varargout{1}=nLicks;
            end
        end
        function zMatrix=rewardHistoryMat(obj,nTrials,bothCS,binary)
            %Create a matrix of reward history where both columns and rows
            %go backward in time
            %INPUT
            %nTrials: # of trials needed to go backwards. it will dictate
            %   the # columns for the output
            %bothCS: use CS+ and CS- in chronological order. if false, uses
            %   only CS+ trials
            %binary: set to true to use only 1s in the trials that have
            %   rewards. Set to false to use the # of rewards in the trials
            %   and windows
            %OUTPUT
            %zMatrix: If you use CS+ trials only:
            %   1st column is for the pre-CS interval. 2nd column is
            %   all trials. 3rd column is trial 2 until the end. 4th column
            %   is trial 3 until the end. And so on
            %
            %for example on trial 10, the corresponding row of zMatrix will
            %have the pre-CS reward of trial 10, then trial 9 reward, trial
            %8 reward, and so on. The last trial never gets used and every
            %CS column will have at least 1 NaN.
            %   If you use all trials, zMatrix is the same except the
            %   columns are pushed back by 1 and the 1st column indicates
            %   whether the trial is a CS+ (1) or CS- (0)
            if binary
                %only use the binary versions
                [z,zBefore]=rewardHistory(obj); %calculate reward history
            else
                %use the reward counts
                [~,~,z,zBefore]=rewardHistory(obj);
            end
            if bothCS
                %initialize a matrix of NaNs. rows need to account for all
                %trials. columns only go as far back as needed
                zMatrix=NaN(length(z.allCS),nTrials+1);
                %flip the order trials so that it goes backwards in time
                z.allCS=flip(z.allCS);
                %assign the pre-CS variable to the 1st column
                zMatrix(:,1)=flip(zBefore.allCS);
                for back=1:nTrials %for each trial back in time
                    col=back+1; %the corresponding column is 1 higher
                    %identify the correct set of rows within the desired columd
                    %of zMatrix
                    zMatrix(1:end-back,col)=z.allCS(1+back:end);
                end
                %create a column vector that identifies the CS+ vs CS-
                [plusT,minusT]=obj.getTrialTimes(1);
                %convert the 2nd column of CS+ to 1s
                plusT(:,2)=ones(size(plusT(:,2)));
                %convert the 2nd column of CS- to 0s
                minusT(:,2)=zeros(size(minusT(:,2)));
                %collate and sort all CS onset and offset times
                allCS=sortrows([plusT;minusT]);
                %push the column of 1s and 0s into the matrix
                zMatrix=[allCS(:,2) zMatrix]; %horizontal concatenation
            else %only the CS+ trials
                %initialize a matrix of NaNs. rows need to account for all
                %trials. columns only go as far back as needed
                zMatrix=NaN(length(z.CSplus),nTrials+1);
                %flip the order trials so that it goes backwards in time
                z.CSplus=flip(z.CSplus);
                %assign the pre-CS variable to the 1st column
                zMatrix(:,1)=flip(zBefore.CSplus);
                for back=1:nTrials %for each trial back in time
                    col=back+1; %the corresponding column is 1 higher
                    %identify the correct set of rows within the desired columd
                    %of zMatrix
                    zMatrix(1:end-back,col)=z.CSplus(1+back:end);
                end
            end
            %at the back of the matrix concatenate satiety
        end
        function [z,zBefore,varargout]=rewardHistory(obj,period)
            %Determine whether reward was consumed for each CS and ITI
            %INPUTS
            %period: optional. duration from CS onset to consider for calculation
            %   It must be a scalar in ms
            %OUTPUT
            %z: struct with 3 fields: CS+, CS-, ITI. Each field is an
            %   indicator variable for reward. 1 if reward on that trials,
            %   otherwise 0
            %zBefore: struct with 2 fields: CS+, allCS. Each field is an
            %   indicator variable for reward in the "period" before the CS
            %OPTIONAL OUTPUT
            %nonbinarized versions of z and zBefore
            if ~exist('period','var') || isempty(period)
                period=8000; %default
            end
            %grab the time stamps for Ensure and CS's
            ensureTimes=obj.getTimes(7);
            [plusT,minusT]=obj.getTrialTimes(period);
            %count the # deliveries during each CS presentation
            nRewards.CSplus =BehaviorData.calcEvents(ensureTimes,plusT(:,1),plusT(:,2));
            nRewards.CSminus =BehaviorData.calcEvents(ensureTimes,minusT(:,1),minusT(:,2));
            %collate and sort all CS onset and offset times
            allCS=sortrows([plusT;minusT]);
            %create a field for all CS collated
            nRewards.allCS=BehaviorData.calcEvents(ensureTimes,allCS(:,1),allCS(:,2));
            %the ITI onsets are the CS offsets
            ITIonset=[0; allCS(:,2)]; %include the time before first CS
            %the ITI offsets are the CS onsets
            ITIoffset=[allCS(:,1); obj.duration]; %include time after last CS
            durations=ITIoffset-ITIonset; %ms
            %identify the intervals that are longer than the CS intervals
            indices=find(durations'>period); %row vector of numeric indices
            %subdivide the ITIs so that they are all the same size
            newIntervals=[];
            for i=indices %for each oversized interval
                %numerate a new row of intervals that counts backward from
                %the end of the ITI
                newRow=ITIoffset(i):-period:ITIonset(i);
                %flip this row so that it moves forward in time
                newRow=flip(newRow);
                %pair up the numbers into onsets and offsets
                newColumns=[newRow(1:end-1)' newRow(2:end)'];
                %append the new columns to the larger matrix
                newIntervals=[newIntervals; newColumns];
            end
            %all ITIs that were shorter than period will be omitted
            nRewards.ITI=BehaviorData.calcEvents(ensureTimes,newIntervals(:,1),newIntervals(:,2));
            %binarize the trials because we only care about the presence of
            %reward not magnitude
            z.CSplus=nRewards.CSplus>0;
            z.CSminus=nRewards.CSminus>0;
            z.allCS=nRewards.allCS>0;
            z.ITI=nRewards.ITI>0;
            %for a conditioning session
            if strcmp(obj.behavior,'conditioning')
                %modify the counting for the CS+ such that there has to be
                %a lick within 1 sec of reward delivery
                %first calculate the latency to lick after the Ensure deliveries
                latencies = getLatencies(obj,7,6,2000,true);
                z.CSplus=latencies<1;
                nRewards.CSplus=int(z.CSplus);
            end
            %initialize the count for pre-CS intervals
            nRewardsBefore.CSplus=zeros(size(z.CSplus));
            for n=1:size(plusT,1) %for each CS+ trial
                %determine the no cue interval that precedes this CS presentation
                %by using the smallest negative difference from this CS onset time
                intervalIndex=find(newIntervals(:,1)-plusT(n,1)<0,1,'last');
                %assign this ITI "trial" into the new struct
                nRewardsBefore.CSplus(n)=nRewards.ITI(intervalIndex);
            end
            nRewardsBefore.allCS=zeros(size(allCS,1),1);
            for n=1:size(allCS,1) %for each trial
                %use the same approach but now without regard for whether
                %it is a CS+ or CS-
                intervalIndex=find(newIntervals(:,1)-allCS(n,1)<0,1,'last');
                nRewardsBefore.allCS(n)=nRewards.ITI(intervalIndex);
            end
            %create binary versions
            zBefore.CSplus=nRewardsBefore.CSplus>0;
            zBefore.allCS=nRewardsBefore.allCS>0;
            if nargout>2
                varargout{1}=nRewards;
                if nargout==4
                    varargout{2}=nRewardsBefore;
                end
            end
        end
        function [totalForCS, noCue]=volumePerCS(obj,period,flag)
            %Calculates the total uL delivered after CS onset for all CS presentations.
            %INPUTS
            %period: duration from CS onset to consider for calculation, in
            %   ms. It must be a scalar or a vector
            %flag: if true, calculates the average instead of the total
            %OUTPUTS
            %totalForCS: 1x2 vector, mL. column1 is CS+, column2 is CS-
            %            if a vector of periods is given as input, there
            %            will be a row for each of those periods
            %noCue: total volume for intervals without cues, mL
            % Created by Francisco Pena on 9/4/18
            
            volumePerDel=volPerD(obj)/1000; %convert from uL to mL
            ensureTimes=obj.getTimes(7);
            [plus,minus]=obj.getTrialTimes(period(1));
            if isscalar(period)
                %count the # deliveries during each CS presentation
                nDeliveriesCSplus =BehaviorData.calcEvents(ensureTimes,plus(:,1),plus(:,2));
                nDeliveriesCSminus =BehaviorData.calcEvents(ensureTimes,minus(:,1),minus(:,2));
                %calculate the sums across trials
                totalForCS=[sum(nDeliveriesCSplus) sum(nDeliveriesCSminus)]*volumePerDel; %mL
                if nargout==2
                    %calculate this before the totalForCS becomes an average
                    noCue=obj.volForTest - sum(totalForCS);
                end
                if nargin==3 && flag
                    %calculate the average instead of the total
                    totalForCS=[mean(nDeliveriesCSplus) mean(nDeliveriesCSminus)];
                    totalForCS=totalForCS*volPerD(obj); %uL
                end
                
            else %for vectors and arrays
                %check that period is a column vector
                if size(period,1)>1 %does it have multiple rows?
                    if size(period,2)==1 %only 1 column
                        period=period'; %convert it to a row
                    else %if it has more than 1 row, it's an array
                        error('input must be a scalar or vector, not an array')
                    end
                end
                nOffsets=length(period);
                %create matrices of offsets: rows are trials, each column is a
                %set of offsets based on a duration
                plusOffsets=repmat(plus(:,1),1,nOffsets)+repmat(period,size(plus,1),1);
                minusOffsets=repmat(minus(:,1),1,nOffsets)+repmat(period,size(minus,1),1);
                %initialize matrices
                nDeliveriesCSplus=zeros(size(plusOffsets));
                nDeliveriesCSminus=zeros(size(minusOffsets));
                for j=1:nOffsets %for each duration, populate the columns
                    %count the # deliveries during each CS presentation
                    nDeliveriesCSplus(:,j) =BehaviorData.calcEvents(ensureTimes,plus(:,1),plusOffsets(:,j));
                    nDeliveriesCSminus(:,j) =BehaviorData.calcEvents(ensureTimes,minus(:,1),minusOffsets(:,j));
                end
                %calculate the sums across trials, which is the 1st dimension
                totalForCS=[sum(nDeliveriesCSplus)' sum(nDeliveriesCSminus)']*volumePerDel;
                if nargin==3 && flag
                    %calculate the average instead of the total
                    totalForCS=[mean(nDeliveriesCSplus)' mean(nDeliveriesCSminus)']*volumePerDel;
                    %convert from mL to uL
                    totalForCS=totalForCS*1000; %uL
                end
            end
            
        end
        function plotVforCS(obj,durations,V)
            %Plots the output of volumePerCS, showing a line for each trial
            
            figure; hold on
            %CS+ in teal
            plot(durations',V(:,1),'LineWidth',2,'Color',[0 0.7 0.6])
            %CS- in yellow
            plot(durations',V(:,2),'LineWidth',2,'Color',[0.9 0.9 0])
            ylim([0 25]); ylabel('cumulative consumption, uL')
            xlabel('time from CS onset')
            xlim([0 max(durations)])
            title(obj.graphTitle)
        end
        function [nLicksCS,CSbool,nLITI,ITIdur]=dataForModel(obj,period)
            %Retrive the data used for fitting a model on a trial-by-trial
            %basis
            %INPUT
            %period:time from CS onset to use for calculation, ms
            %OUTPUT
            %nLicksCS: column vector of lick counts. each row is a trial
            %CSbool: column vector that says whether a trial is CS+ or not
            %nLITI: column vector of lick counts for all ITI
            %ITIdur: column vector listing the duration of all ITIs in sec
            [nLicksCS,~,CSbool]=obj.licksAndDel(period);
            [nLITI,ITIdur]= obj.calcEventsInITI(6,period);
            ITIdur=ITIdur/1000; %convert ms to sec
        end
        function plotLickCounts(obj,period)
            %Shows violin plot for all trials and ITIs
            %INPUT
            %period:time from CS onset to use for calculation, ms
            [nLicksCS,CSbool,nLITI,~]=obj.dataForModel(period);
            %create struct b/c each category may have a different length
            s.CSplus=nLicksCS(CSbool);
            s.CSminus=nLicksCS(not(CSbool));
            s.ITI=nLITI;
            figure
            violinplot(s);
            ylabel('# licks')
            title(obj.graphTitle)
        end
        function [nLicks,nDel,varargout]=licksAndDel(obj,period)
            %Counts licks and deliveries in all CS presentations
            %
            %INPUT
            %period:time from CS onset to use for calculation, ms
            %   it can also be a 2-element vector to indicate start and
            %   stop limits
            %OUTPUTS
            %nLicks: column vector of lick counts. each row is a trial
            %nDel: column vector. # ensure deliveries
            %bool: column vector that says whether a trial is CS+ or not
            %allT: 3 columns: CS onset, CS offset. for all CS
            %   within the function, these times are in ms. the output is sec
            
            lickTimes=obj.getTimes(6); %grab the lick times
            ensureT=obj.getTimes(7); %times of ensure deliveries
            [plus, minus]=obj.getTrialTimes(max(period)); %CS onset times
            %use a new column to indicate whether the trial is a CS+ or not
            plus=[plus ones(size(plus,1),1)];
            minus=[minus zeros(size(minus,1),1)];
            %collect the times and place them in chronological order
            allT=sortrows([plus;minus]);
            if length(period)==2 %if there are 2 entries for period
                %use the 1st one to change the "onset" times
                allT(:,1)=allT(:,1)+period(1);
            end
            %count # licks during each CS presentation
            nLicks=BehaviorData.calcEvents(lickTimes,allT(:,1),allT(:,2));
            %count # ensure deliveries
            nDel=BehaviorData.calcEvents(ensureT,allT(:,1),allT(:,2));
            if nargout==3 %if the user asks for the bool
                varargout{1}=allT(:,3)==1; %create the boolean
            end
            allT=allT(:,1:2)/1000; %convert ms to sec
        end
        function h=plotLickAndDel(obj,period,nLicks,bool,nDel)
            %Plots the #licks and #deliveries as a function of trial #
            %INPUT
            %period: duration from CS onset to use for calculation, in ms
            %OPTIONAL INPUT
            %nLicks: column vector of lick counts. each row is a trial
            %bool: column vector that says whether a trial is CS+ or not
            %nDel: column vector. # ensure deliveries
            %OUTPUT
            %h:handle for the graph
            
            if nargin==2 %if the lick counts were not provided
                %calculate them
                [nLicks,nDel,bool]=obj.licksAndDel(period);
                [nLicksITI,durITI]= obj.calcEventsInITI(6,period);
            end
            nLicksPlus = nLicks(bool);
            nLicksMinus = nLicks(not(bool));
            if ~exist('nDel','var') %if nDel was not calculated or provided
                nDel=zeros(size(nLicks)); %create a filler
            end
            nDelPlus=nDel(bool);
            nDelMinus=nDel(not(bool));
            X=1:size(nLicks,1); %x coordinates
            
            h=figure; %graph for both CS+ and CS-
            %plot the # licks for each trial as a function of trial number
            subplot(2,1,1); hold on
            Y=nLicks/(period/1000); %licks/sec
            %leave CS- in black
            plot(X,Y,'o--k') %dashed line to connect dots
            %CS+ in teal
            scatter(X(bool),Y(bool),'filled','MarkerFaceColor',[0 0.7 0.6])
            ylabel('licks/sec')
            xlabel('trial number')
            legend('CS-','CS+','Location','Best')
            title([obj.graphTitle ' all trials'])
            %ITI plot
            subplot(2,1,2); hold on
            Y=nLicksITI./(durITI/1000); %licks/sec
            %leave CS- in black
            plot(Y,'o--k') %dashed line to connect dots
            ylabel('licks/sec')
            xlabel('ITI number')
            currentTitle=[obj.graphTitle ' trialHistory'];
            savefig(h,currentTitle)
            saveas(gcf,[currentTitle '.png'])
            
            h=figure;
            %CS+ trials on column 1
            subplot(2,2,1); hold on
            plot(nLicksPlus,'o--','color',[0 0.7 0.6])
            ylabel('# licks')
            xlabel('trial number')
            title([obj.graphTitle ' CS+'])
            subplot(2,2,3); hold on
            plot(nDelPlus,'o--','color',[0 0.7 0.6])
            ylabel('# deliveries')
            xlabel('trial number')
            %CS- trials on column 2
            subplot(2,2,2); hold on
            plot(nLicksMinus,'o--k') %dashed line to connect dots
            ylabel('# licks')
            xlabel('trial number')
            title([obj.graphTitle ' CS-'])
            subplot(2,2,4); hold on
            plot(nDelMinus,'o--k') %dashed line to connect dots
            ylabel('# deliveries'); ylim([0 2])
            xlabel('trial number')
            currentTitle=[obj.graphTitle ' trialHistory grouped by CS'];
            savefig(h,currentTitle)
            saveas(gcf,[currentTitle '.png'])
        end
        function R=calcLickRates2(obj,period,nLicksP,nLicksM,nLicksITI)
            %Calculate the mean lick rate for CS+, CS-, and ITI
            %INPUT
            %period: duration from CS onset to use for calculation, in ms
            %nLicksP: matrix where rows are CS+ trials. columns are experiments
            %nLicksM: same for CS-
            %nLicksITI: same for intertrial intervals
            %OUTPUT
            %R: matrix with 3 rows; one for CS+, CS-, ITI. licks/sec.
            %   columns represent different experiments
            totalTimes=[size(nLicksP,1); size(nLicksM,1)]*period/1000; %sec
            %all remaining time is considered ITI time
            ITItime=obj.duration/1000-sum(totalTimes); %sec
            totalTimes=[totalTimes; ITItime]; %column vector, 3x1
            totalTimes=repmat(totalTimes,1,size(nLicksP,2)); %repeat the column for each experiment
            %element-wise division to calculate the average lick rates
            R=[sum(nLicksP); sum(nLicksM); sum(nLicksITI)]./totalTimes;
        end
        function [R, varargout]=calcLickRates(obj,period,nLicks,bool,nLITI)
            %Calculates the ave lick rate after CS onset for all CS presentations.
            %
            %INPUT
            %period: duration from CS onset to use for calculation, in ms
            %OPTIONAL INPUTS
            %nLicks: column vector of lick counts. each row is a trial
            %bool: column vector that says whether a trial is CS+ or not
            %nLITI: same nLicks but for inter-trial intervals but its more
            %       likely to be a matrix with multiple columns
            %OUTPUT
            %R:1x3 vector, for CS+, CS-, ITI. licks/sec. Or multiple
            %       rows if simulating multiple experiments
            %OPTIONAL OUTPUTS
            %delta: similar to R except its the change in lick rate from
            %       the 1st half to the 2nd half. negative values mean the
            %       lick rate decreased in the 2nd half
            %prop: proportion of trials that have at least 1 lick. this
            %      output cannot be called w/o calling for delta
            
            if nargin==2 %if the lick counts were not provided
                %calculate them
                [nLicks,~,bool]=obj.licksAndDel(period);
                [nLITI]= obj.calcEventsInITI(6,period);
            end
            nSim=size(nLicks,2);
            %check that nLicks and bool have the same number of columns
            assert(nSim==size(bool,2),'nLicks and bool do not have the same # columns')
            [R,delta,prop]=deal(zeros(nSim,3)); %initialize variables
            for i=1:nSim %for each simulation
                %column vectors
                nLicksPlus = nLicks(bool(:,i),i);
                nLicksMinus = nLicks(not(bool(:,i)),i);
                nLicksITI=nLITI(:,i); %column
                %all remaining licks are considered to be in the ITI
%                 ITItotal=obj.getNum(6)-sum(nLicksPlus)-sum(nLicksMinus);
%                 miss=ITItotal-sum(nLicksITI); %look at the diff
%                 if miss>0 %if there are licks that were not accounted for
%                     warning('%d licks were left out',miss)
%                 end
                totalTimes=[length(nLicksPlus) length(nLicksMinus)]*period/1000; %sec
                %all remaining time is considered ITI time
                ITItime=obj.duration/1000-sum(totalTimes); %sec
                totalTimes=[totalTimes ITItime]; %row vector, 1x3
                %element-wise division to calculate the average lick rates
                R(i,:)=[sum(nLicksPlus) sum(nLicksMinus) sum(nLicksITI)]./totalTimes;
                if nargout>=2 %if the 2nd output is called
                    %calculate delta for each vector
                    delta(i,1)=BehaviorData.calcDelta(nLicksPlus);
                    delta(i,2)=BehaviorData.calcDelta(nLicksMinus);
                    delta(i,3)=BehaviorData.calcDelta(nLicksITI);
                    %delta calculation would be more efficient if I could
                    %feed a whole matrix of simulations to the function
                    if nargout==3 %if the 3rd output is called
                        %calculate the proportion of engagement
                        prop(i,1)=mean(nLicksPlus>0);
                        prop(i,2)=mean(nLicksMinus>0);
                        prop(i,3)=mean(nLicksITI>0);
                    end
                end
            end
            if nargout>=2
                varargout{1}=delta;
                if nargout==3
                    varargout{2}=prop;
                end
            end
        end
        function inter=showInterLick(obj,nBins)
            %plots the interlick intervals
            %INPUT
            %nBins: # bins for histogram. optional
            %OUTPUT
            %inter: intervals between adjacent licks, ms
            
            if nargin==1 || isempty(nBins) %if argument wasn't provided
                nBins=25; %default value
            end
            inter=diff(obj.getTimes(6)); 
            figure;
            %convert ms to sec
            histogram(inter/1000,nBins,'Normalization','probability')
            xlabel('inter-lick interval, sec')
            ylabel('fraction')
            nLicks=obj.getNum(6);
            title([obj.graphTitle ', nLicks=' num2str(nLicks)])
        end
        function [varargout]=lickHistogram(obj,period,nbins,flag,nLicks,bool)
            %Plots a histogram of #licks on CS+ and CS- trials
            %INPUT
            %period: duration from CS onset to consider, in ms
            %nbins: # bins for fitting the distribution
            %OPTIONAL INPUTS
            %flag: true to graph, false to avoid graphing. default is false
            %nLicks: column vector or column matrix of lick counts. each row is a trial.
            %bool: column vector or column matrix that says whether a trial is CS+ or not
            %OUTPUTS ARE ALL OPTIONAL
            %all 3 outputs are structs with a field for CS+ and a field for
            %CS-
            %lambda: parameter for poisson fit
            %p: probability of rejecting the null hypothesis that the data
            %   comes a poisson distribution
            %chi:struct of stats from the chi squared goodness of fit
            
            if nargin<5 %if the lick counts were not provided
                %calculate them
                [nLicks,~,bool]=obj.licksAndDel(period);
            end
            if nargin<4 || isempty(flag) %if flag was not provided
                flag=false;
            end
            nSim=size(nLicks,2);
            %check that nLicks and bool have the same number of columns
            assert(nSim==size(bool,2),'nLicks and bool do not have the same # columns')
            [lambda,p]=deal(struct('plus',zeros(nSim,1),'minus',zeros(nSim,1)));
            chi=struct('plus',[],'minus',[]); %values will be structs themselves
            for i=1:nSim %for each simulation
                nLicksPlus = nLicks(bool(:,i),i);
                nLicksMinus = nLicks(not(bool(:,i)),i);
                if flag
                    h=figure; hold on
                    %CS+ in teal
                    histogram(nLicksPlus,'BinMethod','integers','FaceColor',[0 0.7 0.6])
                    %CS- in yellow
                    histogram(nLicksMinus,'BinMethod','integers','FaceColor',[0.9 0.9 0])
                    xlabel('# licks')
                    ylabel('# trials')
                    legend('CS+','CS-','Location','NorthEast')
                    title(obj.graphTitle)
                    savefig(h,[obj.graphTitle ' lickNum histo'])
                end
                if nargout>0 %if you ask for outputs
                %test whether the distributions of lick counts are poisson
                [pd,p.plus(i),chi.plus]=BehaviorData.testPoisson(nLicksPlus,nbins);
                try
                    lambda.plus(i)=pd.lambda;
                catch ME 
                    %if the distribution was not poisson
                    lambda.plus(i)=NaN;
                end
                [pd,p.minus(i),chi.minus]=BehaviorData.testPoisson(nLicksMinus,nbins);
                %exponential distribution may be fitted if poisson failed
                switch pd.DistributionName
                    case 'Poisson'
                        lambda.minus(i)=pd.lambda;
                    case 'Exponential'
                        lambda.minus(i)=pd.mu;
                end
                varargout{1}=lambda;
                    if nargout>1
                        varargout{2}=p;
                        if nargout>2
                            varargout{3}=chi;
                        end
                    end
                end
            end
            
        end        
        function output=calcAnticipation(obj,period,threshold)
            %Calculates the fraction of trials with anticipatory licks. A
            %trial only counts as anticipatory if at least x licks occur,
            %where x=threshold.
            %OPTIONAL INPUT
            %period: duration from CS onset to consider for calculation, in
            %   ms.
            %threshold: # licks required for a trial to count as
            %   anticipatory
            %OUTOUT
            %output:1x2 vector. 1st entry is for the CS+, 2nd is CS-
            % Created by Francisco Pena on 9/14/18
            if nargin<3 || isempty(period) %if period is not provided
                %default is to include the cue and the trace interval
                period=obj.cueDurationMS+obj.traceMS;
            end
            %if threshold was not provided
            if ~exist('threshold','var') || isempty(threshold)
                threshold=1; %default
            end
            p=NaN; m=NaN;
            %for CS-, the interval for anticipation includes the CS and the
            %mean trace interval
            [plus,minus]=obj.getTrialTimes(period);
            nDeliveries=obj.getNum(7); %works for chocolate or vanilla
            lickTimes=obj.getTimes(6); %grab the lick times
            if ~isempty(plus) %if CS+ presentations occurred
                %for CS+, the trace period is variable so use the Ensure times to
                %indicate the end of the trace interval
                if size(plus,1) ==nDeliveries %if there's an equal # trials and deliveries
                    plus(:,2)=obj.getTimes(7); %works for chocolate or vanilla
                elseif size(plus,1) == (nDeliveries-1) %there's 1 CS+ presentation w/o a delivery
                    plus(1:end-1,2)=obj.getTimes(7);
                %if there were no rewards nad only CS+ presentations,
                %use the offset time already calculated above
                end
                %count the # licks during each CS presentation
                nLicksPlus =BehaviorData.calcEvents(lickTimes,plus(:,1),plus(:,2));
                %calculate the proportion of trials w/a lick during the
                %anticipation period
                p=sum(nLicksPlus>=threshold)/length(plus);
            end
            if ~isempty(minus) %if CS- presentations occurred
                %count the # licks during each CS presentation
                nLicksMinus =BehaviorData.calcEvents(lickTimes,minus(:,1),minus(:,2));
                m=sum(nLicksMinus>=threshold)/length(minus);
            end
            output=[p m];
        end
        function [counts,varargout] = calcEventsInITI(obj,label,period)
            %Calculate the #events during intertrial intervals
            %INPUTS
            %label: an integer that refers to a label for timestamps
            %period: in ms
            %OUTPUT
            %counts: vector with the # for each ITI
            %dur: durations of each interval in ms
            eventTimes=obj.getTimes(label);
            [plus, minus]=obj.getTrialTimes(period);
            %collate and sort all CS onset and offset times
            allCS=sortrows([plus;minus]);
            %the ITI onsets are the CS offsets
            ITIonset=[0; allCS(:,2)]; %include the time before first CS 
            %the ITI offsets are the CS onsets
            ITIoffset=[allCS(:,1); obj.duration]; %include time after last CS
            counts=BehaviorData.calcEvents(eventTimes,ITIonset,ITIoffset);
            if nargout==2 
                %calculate durations
                varargout{1}=ITIoffset-ITIonset;
            end
        end
        function N = calcEventsBeforeCS(obj,label,dur)
            %INPUTS
            %label: 6, 7, or 8 for licks or Ensure deliveries
            %dur: duration in ms to use before the CS
            %OUTPUT
            %N: struct with a field for CS+ that is a vector of # events
            %   for each pre-trial interval. There will be an analagous
            %   field for the CS- if there are CS- presentations.
            
            eventTimes=obj.getTimes(label);
            %collect the CS onset times
            [plus, minus]=obj.getTrialTimes(-dur);
            if ~isempty(plus) %if there are CS+ trials
                %the second column has the timestamps that begin before the CS
            N.plus=BehaviorData.calcEvents(eventTimes,plus(:,2),plus(:,1));
%             fprintf('%d of %d CS+ trials have licks\n',sum(N.plus>0), numel(N.plus))
            end
            if ~isempty(minus) %if there are CS- trials
            N.minus=BehaviorData.calcEvents(eventTimes,minus(:,2),minus(:,1));
%             fprintf('%d of %d CS- trials have licks\n',sum(N.minus>0),numel(N.minus))
            end
        end
        %% preparation before constructing ImagingAndBehavior objects
        function [nMissing,onsets,obj]= checkForMissing(obj)
            %Check if there are any timestamps for Ensure missing. The
            %assumption only works for session of conditioning, not
            %consumption tests. However, a different algorithm could be
            %applied to those tests by using the FR schedule.
            %All times are in ms
            %OUTPUTS
            %nMissing: # stamps that were missing
            %onsets: times of CS+ onset that do not have a matching delivery
            %obj:updated object
            switch obj.behavior
                case {'CS+ only', 'conditioning','Conditioning'}
                    %acceptable cases
                otherwise
                    error(['data is not from a conditioning session. behavior: ' obj.behavior])
            end
            if obj.getNum(-8)>0 %if this data has already been checked and given a warning stamp
                fprintf('There are %d warning stamps already in the file',obj.getNum(-8))
            end
            nMissing=0; onsets=[]; %initialize the outputs
            plusTimes=obj.getTimes(obj.getlabelOfCS('CS+'));
            nPlus=length(plusTimes);
            nDel=obj.getNum(7); %no. deliveries
            %what's the maximum time between CS onset and a delivery?
            maxDur=6000; %ms
            %for each CS+ presentation there should be an Ensure delivery
            if nPlus>nDel %if there's more CS's
                ensTimes=obj.getTimes(7);
                %determine which CS+ presentations do not have a matching Ensure
                for i=1:(nPlus-1) %for each cue, except the last one
                    %sometimes the last delivery gets cut off because time runs out
                    onset=plusTimes(i);
                    allDiff=ensTimes-onset; %times relative to 1 cue
                    %is there a duration within the expected window?
                    if any(allDiff(allDiff>0)<maxDur) %positive and less than the max
                        %great
                    else %if not
                        nMissing=nMissing+1; %increment the counter
                        onsets=[onsets;onset]; %add this time to the output
                        %insert a warning stamp for this trial
                        %and insert an approximate timestamp for the delivery
                        approx=onset+obj.cueDurationMS+obj.traceMS;
                        %sort the rows based on 2nd column (time)
                        obj.box=sortrows([obj.box;-8 onset; obj.labelOfFlavor approx],2);
                    end
                end
                %the output object will have the updated 'box' property
                saveObj(obj);
            elseif nPlus <nDel %if there's less
                %not sure why this would happen in a conditioning session
                warning(['There are more deliveries than CS+ presentations for ' obj.graphTitle])
            end %otherwise, they're equal so there's no problem
        end
        function [putBack, other]=parseInterrupted(obj)
            % Extracts the easy cases of the interrupted data into 2 outputs
            %   obj.box:struct with a field called interruptedData
            %   interruptedData: cell array
            %   putBack: 2 column matrix of labels and timestamps. ready to put into
            %           the matrix of regular data
            %   other: cell array of data that did not fit an easy case. It needs to be
            %           parsed on a case-by-case basis
            
            if isstruct(obj.box)
                interruptedData=obj.box.interruptedData;
            else
                error('There is no interrupted data for this file')
            end
            possibleLabels=[1 2 6 7 8];
            putBack=[];
            other={};
            %only proceed if there is data to be parsed
            if ~isempty(interruptedData)
                for i=1:length(interruptedData) %for each row
                    currentRow=cell2mat(interruptedData(i));
                    %rows with only 1 element don't offer enough info
                    if length(currentRow) ==3 %only analyze the rows that have 3 elements
                        %compare it to an expected 2-label case
                        if any(currentRow(1)==possibleLabels) %check the 1st element
                            if currentRow(2)==10 %check the 2nd element
                                %if it is a 2-label case, check that the next element is a
                                %timestamp? Not sure how I would do that
                                new=[currentRow(1) currentRow(3);currentRow(2:3)];
                                putBack=[putBack;new]; %store the newly arranged data
                            else %if the 2nd label is not a 10 (SYNC pulse)
                                other=[other;currentRow]; %store this row
                                %and the next one if it is a timestamp
                                %sometimes the next row is a row with 3 elements
                                if i+1 <=length(interruptedData) %check that the next row exists
                                    nextRow=cell2mat(interruptedData(i+1));
                                    if length(nextRow)==1
                                        other=[other;interruptedData(i+1)];
                                    end
                                end
                            end
                        else %if the first element doesn't fit an expected label
                            other=[other;currentRow]; %store this row
                        end
                    end
                end
            else %Otherwise, throw a warning message
                warning('There is no interrupted data for this file')
            end
        end
        %% Analyzing optogenetics data
        function [yes,no]=optoOrNot(obj,whichCS,skips)
            %Identifies which trials had optogenetic stimulation
            %INPUTS
            %whichCS: 'CS+' or 'CS-'
            %skips: a vector of indices for trials that should be skipped
            %   from this analysis. optional
            %OUTPUTS
            %column vectors of boolean indices
            
%             %check for timestamps of the optogenetic stim
%             if obj.getNum(12)<1 %if there are none
%             warning('There was no optogenetic stimulation.')
%             end
            %collect the timestamps
            optoT=obj.getTimes(12);
            CStimes=obj.getTimes(obj.getlabelOfCS(whichCS));
            assert(~isempty(CStimes),'There are no presentations of the desired CS')
            %check which CS times have a corresponding opto time
            numCue=length(CStimes);
            yes=false(numCue,1); %initialize a boolean
            for i=1:numCue %for each cue
                allDiff=optoT-CStimes(i); %times relative to 1 cue
                %if there is a difference less than 2ms then this cue has
                %optogenetic stim
                yes(i)=any(abs(allDiff)<2);
            end
            
            if exist('skips','var') && ~isempty(skips) %if skips are provided
                %remove trials that need to be skipped
                yes(skips)=[];
            end
            no=not(yes);
        end
        function [R,Ropto]=lickRates(obj)
            %Calculates lick rates for all trials and separates those with
            %opto
            %INPUT
            %
            %OUTPUT
            %R:1x3 vector, for CS+, CS-, ITI. average rate, licks/sec.
            %Ropto: 1x2 vector. CS+ and CS- rates during opto trials
            R=zeros(1,3); Ropto=zeros(1,2);
            [~,ratesITI]=createTrialTable(obj);
            R(3)=mean(ratesITI);
            T=trials2summ(obj); %create summary table 
            R(1)=T.meanLickRate(2); %CS+
            R(2)=T.meanLickRate(4); %CS-
            Ropto(1)=T.meanLickRate(1); %CS+ w/opto
            Ropto(2)=T.meanLickRate(3); %CS- w/opto
        end
        function output=forOptoTable(obj)
            %This method should only be used for consumption tests w/o cues
            %create summary table with the following columns
            vNames={'isCS+','hasOpto','meanV','totalV','meanLickRate','meanLatency'};
            if obj.getNum(12)>0 %if there's opto
                %limit the volume to the first 20 min
                timeLimit=20*60; %sec
                V=obj.volDuring(0,timeLimit); %inputs in s
                lickRate=obj.lickDuring(0,timeLimit); %lick/s
            else %otherwise use the full amount
                V=obj.volForTest; %mL
                lickRate=obj.getNum(6)/(obj.duration/1000); %licks/s
            end
            output=table(false,any(obj.getNum(12)),NaN,V,...
                lickRate, NaN,'VariableNames',vNames);
        end
        function [T,rateITI]=createTrialTable(obj,period)
            %Create a table with an entry for every trial
            %OPTIONAL INPUT
            %period: time from CS onset to use for calculations, in ms
            %OUTPUT
            %T: the columns of the table are listed in vNames below
            %rateITI: optional. lick rate during each ITI. column vector
            if ~exist('period','var') || isempty(period)
                %default value
                period=8000; %8 sec in ms
            end
            vNames={'isCS+','hasLick','hasOpto','volume','lickRate','lickLatency','lickCount'}; %names of the columns
            if obj.getNum(1)>0 && obj.getNum(2)>0 %if there are CS presentations
                volumePerDelivery=volPerD(obj); %uL
                %count the # licks and # deliveries during each CS presentation
                [nLicksCS,nDelCS,isCSplus]=obj.licksAndDel(period);
                %isCSplus is a boolean that indexes into the other 2
                %outputs of licksAndDel
                %convert from deliveries to volume
                vCSplus=nDelCS(isCSplus)*volumePerDelivery;
                vCSminus=nDelCS(~isCSplus)*volumePerDelivery;
                %convert lick counts to lick rates
                lickRateCS=nLicksCS/(period/1000); %licks/sec
                if nargout==2
                %count the licks in the ITI as well
                [nLicksITI,dur]= obj.calcEventsInITI(6,period);
                %convert dur from ms to sec
                rateITI=nLicksITI./(dur/1000); %durations of ITI vary
                end
                [yL,~]=obj.lickOrNot('CS+', period);
                [yO,~]=obj.optoOrNot('CS+');
                %calculate the latency to lick relative to CS+
                label=obj.getlabelOfCS('CS+');
                latencyCap=10000; %10 sec in ms
                latencies = getLatencies(obj,label,6,latencyCap,true); %sec
                %put together the CS+ data
                T=table(true(size(yL)),yL,yO,vCSplus,lickRateCS(isCSplus),...
                    latencies,nLicksCS(isCSplus), 'VariableNames',vNames);
                
                [yL,~]=obj.lickOrNot('CS-', period);
                [yO,~]=obj.optoOrNot('CS-');
                label=obj.getlabelOfCS('CS-');
                latencies = getLatencies(obj,label,6,latencyCap,true); %sec
                %put together the CS- data
                temp=table(false(size(yL)),yL,yO,vCSminus,lickRateCS(~isCSplus),...
                    latencies,nLicksCS(~isCSplus),'VariableNames',vNames);
                T=[T;temp]; %concatenate the tables
            else %if there are no CS presentations
                %provide a default table
                T=table(false,false,false,NaN,NaN,NaN,'VariableNames',vNames);
            end
        end
        function output=trials2summ(obj,period)
            %Summarize the optogenetics trials into 4 behavior measures
            %OPTIONAL INPUT
            %period: time from CS onset to use for calculations, in ms
            %OUTPUT
            %output: a table that summarizes the trials based on CS and
            %   whether there is opto. the last 2 columns are: 
            %lickprop: proportion of trials with a lick
            %latency: mean lick latency for trials with licks
            if ~exist('period','var')
                period=[]; %default
            end
            %create a table for all trials
            trialT=createTrialTable(obj,period);
            period=period/1000; %convert period from ms to s
            %create summary table with the following columns
            vNames={'isCS+','hasOpto','meanV','totalV','meanLickRate',...
                'meanLatency','lickprop','latency'}; 
            nRows=4;
            blank=zeros(nRows,1);
            %initialize table and populate the first 2 columns
            output=table(logical([1 1 0 0]'),logical([1 0 1 0]'),...
                blank,blank,blank,blank,blank,blank,'VariableNames',vNames);
            %create boolean columns that index into the trial table
            trialType=cell(1,4); %intialize a cell array to keep them in
            trialType{1}=trialT.('isCS+') & trialT.hasOpto;
            trialType{2}=trialT.('isCS+') & not(trialT.hasOpto);
            trialType{3}= ~trialT.('isCS+') & trialT.hasOpto;
            trialType{4}=~trialT.('isCS+') & ~trialT.hasOpto;
            %populate the blank columns of the new table
            for i=1:nRows %for each row of new table
                %calculate the mean volume 
                output{i,3}=mean(trialT.volume(trialType{i}));
                %calculate the total volume
                output{i,4}=sum(trialT.volume(trialType{i}));
                %calculate the mean lick rate
                output{i,5}=mean(trialT.lickRate(trialType{i}));
                %query for the latencies of this trial type
                lat=trialT.lickLatency(trialType{i});
                %calculate the mean of all trials
                output{i,6}=mean(lat);
                %proportion of trials that have a lick
                output{i,7}=mean(lat<period);
                if output{i,7}==0 %if there are no licks in the interval
                    output{i,8}=period; %use the maximum instead of a NaN
                else
                    %calculate the mean of trials that have a lick
                    output{i,8}=mean(lat(lat<period));
                end
            end
        end
        function [V, Vopto, noCue]=volumePerCSopto(obj,period,flag)
            %Calculates the total or mean uL delivered after CS onset for all CS presentations.
            %INPUTS
            %period: duration from CS onset to consider for calculation, in
            %   ms. It must be a scalar or a vector
            %flag: if true, calculates the average instead of the total
            %OUTPUTS
            %V: 1x2 vector, column1 is CS+, column2 is CS-
            %       if a vector of periods is given as input, there will be
            %       a row for each of those periods
            %Vopto: same as V but for the trials that have optogenetic stim
            %noCue: total volume for intervals without cues
            %NOTE: all volumes are in uL, not mL
            V=[]; Vopto=V; %initialize empty matrices
            if isscalar(period)
                trialT=obj.createTrialTable(period); %collect and sort data by trial type
                %create boolean columns that index into the trial table
                trialType=cell(1,4); %intialize a cell array to keep them in
                trialType{1}=trialT.('isCS+') & trialT.hasOpto;
                trialType{2}=trialT.('isCS+') & not(trialT.hasOpto);
                trialType{3}= ~trialT.('isCS+') & trialT.hasOpto;
                trialType{4}=~trialT.('isCS+') & ~trialT.hasOpto;
                V(1)=sum(trialT.volume(trialType{2})); %CS+
                V(2)=sum(trialT.volume(trialType{4})); %CS-
                Vopto(1)=sum(trialT.volume(trialType{1})); %CS+
                Vopto(2)=sum(trialT.volume(trialType{3})); %CS-
                noCue=obj.volForTest*1000-sum(V)-sum(Vopto);
                if exist('flag','var') && flag
                    %use the mean instead of the sum
                    V=mean(trialT.volume(trialType{2}));
                    V=[V mean(trialT.volume(trialType{4}))];
                    Vopto=mean(trialT.volume(trialType{1}));
                    Vopto=[Vopto mean(trialT.volume(trialType{3}))];
                    %leave noCue the same
                end
            else %for vectors and arrays
                %check that period is a row vector or change it to be so
                if size(period,1)>1 %does it have multiple rows?
                    if size(period,2)==1 %only 1 column
                        period=period'; %convert it to a row
                    else %if it has more than 1 row, it's an array
                        error('input must be a scalar or vector, not an array')
                    end
                end
                for j=period %for each period
                    %only keep the last iteration of noCue
                    [out1, out2, noCue]=volumePerCSopto(obj,j,flag);
                    %vertically concatenate the volumes
                    V=[V;out1];
                    Vopto=[Vopto;out2];
                end
            end
        end
        function plotLickPSTHopto(obj,binWidth)
            %Plots the PSTH of CS+ and CS-, with and without optogenetics
            %INPUTS
            %binWidth: width for sliding window in sec
            if ~exist('binWidth','var') || isempty(binWidth)
                binWidth=1; %default is 1 s
            end
            %calculate the PSTH for each CS
            [psthPlus,optoPlus,t]=obj.lickPSTHopto('CS+',binWidth);
            [psthMinus,optoMinus,~]=obj.lickPSTHopto('CS-',binWidth);
            %collect the whole set of trials, not just the average
            
            %plot against time
            figure
            yMax=1.5;
            grayColor=[0.5 0.5 0.55];
            subplot(1,2,1); hold on
            %CS trials without opto in black
            plot(t,psthPlus,'k','LineWidth',2)
            %CS trials with opto in gray
            plot(t,optoPlus,'Color',grayColor,'LineWidth',2)
            %draw a rectangle for the duration of the CS
            xRect=[0 obj.cueDurationMS obj.cueDurationMS 0]/1000; %sec
            yRect=[0 0 yMax+0.1 yMax+0.1];
            %create a handle for the rectangle
            hRect=fill(xRect,yRect,[0 0.7 0.6]); %CS+ in teal
            hRect.FaceAlpha=0.5; %set the transparency to 50%
            hRect.EdgeColor='none'; %remove the outline
            %small bar for opto LED period
            plot([0 5],[0.01 0.01],'Color',grayColor,'LineWidth',3)
            legend('control','inactivation')
            xlabel('time, s')
            ylabel('licks/s')
            axis tight
            ylim([0 yMax])
            
            subplot(1,2,2); hold on
            %CS trials without opto in black
            plot(t,psthMinus,'k','LineWidth',2)
            %CS trials with opto in gray
            plot(t,optoMinus,'Color',grayColor,'LineWidth',2)
            %draw a rectangle for the duration of the CS
            xRect=[0 obj.cueDurationMS obj.cueDurationMS 0]/1000; %sec
            yRect=[0 0 yMax+0.1 yMax+0.1];
            %create a handle for the rectangle
            hRect=fill(xRect,yRect,[0.9 0.9 0]); %CS- in yellow
            hRect.FaceAlpha=0.5; %set the transparency to 50%
            hRect.EdgeColor='none'; %remove the outline
            %small bar for opto LED period
            plot([0 5],[0.01 0.01],'Color',grayColor,'LineWidth',3)
            legend('control','inactivation')
            xlabel('time, s')
            ylabel('licks/s')
            axis tight
            ylim([0 yMax])
            
            suptitle(obj.graphTitle)
        end
        function [PSTH,optoPSTH,t]=lickPSTHopto(obj,whichCS,binWidth)
            %Compute a lick PSTH with a sliding window
            %INPUTS
            %whichCS: 'CS+' or 'CS-'
            %binWidth: width of the moving window in sec
            %OUTPUTS
            %PSTH: row vector. the unit is licks/s. only for light OFF
            %   trials
            %optoPSTH: only for light ON trials
            %t: vector of time points for each value in PSTH. s, not ms
            
            binsPerSec=10; %10 bins/s so that the step size is 100ms
            %indicate the # seconds before and after CS onset or leave
            %those inputs empty to use the default [-2 10]
            [~,lickTicks,~,t]=obj.genRaster(whichCS,6,12);
            %convert t from ms to s
            t=t/1000;
            %add the durations before and after CS onset to get the total
            numSec=abs(t(1))+t(end);
            %lickTicks has a time bin at each ms. group the ticks into 
            %non-overlapping bins of 100ms
            lickTicks=binMat(lickTicks,binsPerSec*numSec); % 10 bins/s
            %downsample the time vector as well, 100 ms increments
            t=t(1:100:end-1); %omit the last bin b/c binMat often cuts off the last time bin
            %split up the trials based on optogenetics
            [yes,no]=optoOrNot(obj,whichCS);
            %average across trials using the boolean indices
            PSTH=mean(lickTicks(no,:)); %average across rows (along columns)
            optoPSTH=mean(lickTicks(yes,:));
            %take the moving average so that each value is the average lick
            %rate for the set of bins centered at that bin
            nBins=binsPerSec*binWidth;
            PSTH=movmean(PSTH,nBins)*1000;
            optoPSTH=movmean(optoPSTH,nBins)*1000;
            %the values are multiplied by 1000 at the end to bring the
            %units to licks/s. I have used the calcLickRates() method to
            %verify that these rates make sense empirically.
        end
        %% convenient but not that useful
        function string=printTrialTime(obj,whichCS,trialNum)
            %prints the time in terms of minutes and seconds
            %INPUTS
            %trialNum: integer or vector of integers
            %whichCS: 'CS+' or 'CS-'
            onsets=obj.getTimes(obj.getlabelOfCS(whichCS)); %in ms
            for j=1:length(trialNum) %for each indicated trial
            string=sec2hms(onsets(trialNum(j))/1000) %input in sec
            end
        end
        function prop=getCSproportion(obj)
            %Calculate the proportion of trials that are CS+
            [p,m]=obj.getTrialTimes(0);
            nPlus=size(p,1);
            prop=nPlus/(nPlus+size(m,1));
        end
        function t=getDurations(obj,period)
            %Calculate the total time for CS+ trials, CS- trials, and ITIs
            %INPUT
            %period: sec, amount of time after CS onset
            %OUTPUT
            %t: sec. row vector, 3 entires. CS+, CS-, ITI
            [plus, minus]=getTrialTimes(obj,period);
            t=[length(plus) length(minus)]*period; %sec
            t=[t obj.duration/1000-sum(t)];
        end
        function d=rewardDensity(obj)
            %Calculate the reward density for this session
            %OUTPUT:
            %d: rewards per minute
            %NOTE: the reasont this is not a get method is because there
            %may be a reason to use an input e.g. density in a specific
            %window of time
            
            %divide the # of rewards by amount of time in s
            d=obj.getNum(7)/(obj.duration/1000);
            %convert from rewards/s to rewards/min
            d=d*60;
        end
        function varargout=rollingCount(obj,label,window)
            %Computes the rolling count of rewards (or licks)
            %INPUT
            %label: 7 or 8 for rewards. 6 for licks.
            %window: window of time in sec to use for the rolling interval
            %OUTPUT
            %rolling:
            %time: optional. vector of time in sec with an element for
            %   every element in count. They could be plotted together
            
            %collect the timestamps and convert to centiseconds
            stamps=round(obj.getTimes(label)/10); %cs
            %create a vector of time with an element for every centisecond
            buffer=800; %add 8 sec to the front and the back
            t=(stamps(1)-buffer):(stamps(end)+buffer);
            %initialize the indicator variable
            ticks=zeros(size(t));
            %place a 1 at every stamp
            ticks(ismember(t,stamps))=1;
            nBins=window*100; %convert the window to centiseconds
            rolling=movmean(ticks,nBins)*nBins; %multiply so that the final units are # events in window
            varargout{1}=rolling;
            if nargout==2
                varargout{2}=t/100; %convert cs to s
            end
        end
        function v=cumulativeCon(obj)
            %Calculate the cumulative consumption at all time points
            %get the timestamps of Ensure
            rewardTimes=obj.getTimes(7); %ms
            %create a time vector for the whole session
            t=0:(obj.duration+2000); %give 2-seconds extra in case there were extra imaging frames
            v=zeros(size(t)); %initialize
            v(rewardTimes)=volPerD(obj); %uL
            v=cumsum(v)/1000; %convert uL to mL
        end
    end
    methods (Static) %these methods do not take the obj as an input
        function [nEventsInt,groupedEventTimes]=calcEvents(eventTimes,onset,offset)
            %Calculates the # events during an interval. All times are in ms.
            %This code can also be used on licks and lever presses.
            %INPUTS
            %eventTimes:vector of timestamps for event of interest
            %onset:times of interval onsets
            %offset:times of interval offsets
            %OUTPUT
            %nEventsInt:vector with an element for each interval
            %OPTIONAL OUTPUT
            %groupedEventTimes: cell array that has the corresponding event times
            %NOTE: all inputs need to have the same units (e.g. all in ms
            %or all in sec)
            
            nIntervals=length(onset);
            %check that the lengths of the onset and offset vectors are equal
            assert(nIntervals==length(offset),'The lengths of onset and offset are not equal')
            nEventsInt=zeros(nIntervals,1);
            groupedEventTimes=cell(size(nEventsInt));
            for i=1:nIntervals % for each interval
                %find events (e.g. licks) with times after the CS onset
                after=find(eventTimes>=onset(i));
                %find events with times before the CS offset
                before=find(eventTimes<=offset(i));
                %find the events that fit both criteria
                %true where the elements of A are in B and false otherwise
                between=ismember(after, before); %boolean that can index into after
                nEventsInt(i)=sum(between);
                if nargout==2
                eventInd=after(between); %indices that correspond to the eventTimes
                groupedEventTimes{i}=eventTimes(eventInd);
                end
            end
        end        
        function [score]=calcDelta(X)
            %Calculates the difference between the means of the first and
            %second halves of a column
            %INPUT
            %X: column vector
            %OUTPUTS
            %score: scalar or row vector with an entry for each column of X
            %NOTE: this function is capable of handling a matrix input
            
            %input vector/matrix cannot be empty
            assert(~isempty(X),'input is empty')
            nTrials=size(X,1); %columns are expts
            if nTrials==1 %if a scalar was given
                score=0; %provide an output anyway
            else
                midTrial=floor(nTrials/2); %the midpoint of the trials
                %first,second, and score will be scalars or row vectors together
                first=mean(X(1:midTrial,:)); %average across trials
                second=mean(X((midTrial+1):end,:)); % across trials
                score=second-first;
            end
        end
        function [pd,p,st]=testPoisson(input,nbins)
             %Test whether a distribution of counts (e.g. #licks) is poisson
             %INPUT
             %input: a column vector
             %nbins: #bins for a histogram
             %OUTPUTS
             %pd: object for fitted Possion distribution
             %p: p value of chi squared goodness of fit
             %st:struct from chi squared
             
             %check that the input is not a matrix
             assert(size(input,2)==1,'input is a row vector or matrix')
             %bin the data
             [obsCounts,edges]=histcounts(input,nbins,'BinMethod','integers');
             %convert the edges into bin centers
             bins=edges(1:end-1)+diff(edges)/2;
             %Fit a Poisson probability distribution object to the data
%              pd = fitdist(input,'Poisson'); %can't determine the bins this way
             pd = fitdist(bins','Poisson','Frequency',obsCounts');
             %compute the expected count for each bin
             expCounts = sum(obsCounts)*pdf(pd,bins); %the pdf must be evaluated at integers
             %test the null hypothesis that the observed data comes from
             %the fitted distribution
             try
                 [~,p,st] = chi2gof(bins,'Ctrs',bins,...
                     'Frequency',obsCounts, ...
                     'Expected',expCounts,...
                     'NParams',1);
                 if isnan(p) %if there were not enough degrees of freedom
                     %test for an exponential dist
                     pd = fitdist(bins','Exponential','Frequency',obsCounts');
                     expCounts = sum(obsCounts)*pdf(pd,bins);
                     [~,p,st] = chi2gof(bins,'Ctrs',bins,...
                         'Frequency',obsCounts, ...
                         'Expected',expCounts,...
                         'NParams',1);
                 end
             catch ME %if there's an error (potentially dealing with bins)
                 %provide outputs
                 p=NaN; st=ME;
             end
        end
    end
    
end