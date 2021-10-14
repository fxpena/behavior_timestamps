function output = loadBehaviorData(namesOrIDs,sessionNum,hasN,ca,fear)
%Loads matlab files of behavior data
%   INPUTS
%   namesOrIDs: cell array of strings or vector of integers representing ID
%           numbers
%   sessionNum: integer or vector of integers
%   OPTIONAL INPUTS
%   hasN: if true, will load file names that have an N e.g. 'bN24 CSPF5'
%   instead of 'b24 CSPF5'
%   ca: if true, will load ImagingAndBehavior object instead 
%   fear: if true, loads fearRecall data
%   OUTPUT: cell array of behavior data objects
%NOTE: This function assumes the mice are named with a 'b' such as 'b122'.
%   If the mouse/file name has a different letter such as 'N24' then you
%   need to specify this name i.e. 'bN24 CSPF5'
if ~exist('hasN','var') || isempty(hasN)
    hasN=false; %default
end
if ~exist('ca','var') || isempty(ca)
    ca=false; %default
end
if ~exist('fear','var') || isempty(fear)
    fear=false; %default
end

nFiles=length(namesOrIDs);
output=cell(nFiles,1); %initialize the output
for ind=1:nFiles %for each mouse (or each element in the first input)
    if iscell(namesOrIDs) %if its a cell array
        name=namesOrIDs{ind}; %grab the file name
    else %if its a vector, create the file name based on the inputs
        assert(~isempty(sessionNum),'sessionNum was expected as an input')
        sessionNum=round(sessionNum); %ensure that they're integers
        if length(sessionNum)==1
            %create the file name(s)
            if ~hasN %if it doesn't need an N for the name
                name=sprintf('b%d CSPF%d',namesOrIDs(ind),sessionNum);
            else %if an N is desired
                name=sprintf('bN%d CSPF%d',namesOrIDs(ind),sessionNum);
            end
            if fear %load fearRecall data instead
                if ~hasN %if it doesn't need an N for the name
                    name=sprintf('b%d fearRecall%d',namesOrIDs(ind),sessionNum);
                else %if an N is desired
                    name=sprintf('bN%d fearRecall%d',namesOrIDs(ind),sessionNum);
                end
            end
        else %if the 2nd input is a vector too
            %check whether all of the IDs are the same
            assert(length(namesOrIDs)==length(sessionNum),'the inputs were expected to be of equal length')
            if ~hasN %if it doesn't need an N for the name
                name=sprintf('b%d CSPF%d',namesOrIDs(ind),sessionNum(ind));
            else %if an N is desired
                name=sprintf('bN%d CSPF%d',namesOrIDs(ind),sessionNum(ind));
            end
        end
        if ca %load the calcium imaging file
            name=[name ' with calciumTraces'];
        end
    end
    
    try %load the file
        load([name '.mat']); %append the file ending
        % store each file in the cell array for easy access
        if exist('obj','var') %check for the object
            output(ind) = {obj};
        elseif exist('Ca','var') %check for the imaging object
            output(ind)={Ca};
        end
    catch ME
        %throw the most common error. However, the error might be something
        %else
        warning(['the file ' name '.mat probably does not exist in the directory'])
    end
end
disp('Data was successfully loaded.')
end

