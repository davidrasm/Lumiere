function [ output_args ] = convertPhyForBDMM(file_name)
%   Convert a phy tree simulated in Matlab to BDDMM format for BEAST XML

if ~exist('file_name','var')
  file_name = 'fitBD_neutralFitModel_5sites';
end

% Input .mat tree file
phyFileName = strcat(file_name, '_phy.mat');
load(phyFileName);

%beastIndexing = true;

% Output newick tree file
filePrefix = file_name;
newickFileName = strcat(filePrefix, '_newick.tre');
[phy] = writeNewick(phy, newickFileName);

%%%% Write files for other XML input %%%
[phy,tipNames]=getTipNames(phy);

% Output for locations
fileName = strcat(filePrefix, '_bdmmTipStates.txt');
tipStatesFileID = fopen(fileName, 'w');

% Output for tips
fileName = strcat(filePrefix, '_bdmmTipTimes.txt');
tipTimesFileID = fopen(fileName, 'w');

% Output for dummy seqs
fileName = strcat(filePrefix, '_bdmmAlignment.txt');
alignFileID = fopen(fileName, 'w');


tipDates = zeros(length(tipNames),1);
for n = 1:length(tipNames)

    header = tipNames{n};

    % Write sampling state to file
    sampleLoc = phy.nodes{n}.lineStates(1); % number of mutations
    stateString = strcat(header, '=', num2str(sampleLoc), ',');
    fprintf(tipStatesFileID, '%s', stateString);
    fprintf(tipStatesFileID, '%s\r\n', ''); %puts line break after coordinates

    % Write sampling time to file
    splits = regexp(header, '_', 'split');
    date = splits{end};
    tipDates(n) = str2double(date);
    timeString = strcat(header, '=',date,',');
    fprintf(tipTimesFileID, '%s', timeString);
    fprintf(tipTimesFileID, '%s\r\n', '');

    % Write dummy seqs to file for XML
    seq = num2str(phy.nodes{n}.lineSeqs(1,:)); % convert seq to string
    seq = seq(find(~isspace(seq))); % remove whitespace 
    alignString = strcat('<sequence id="seq_', header, '" taxon="', header, '" totalcount="2" value="', seq, '"/>');
    fprintf(alignFileID, '%s', alignString);
    fprintf(alignFileID, '%s\r\n', '');
        
end

fclose(tipStatesFileID);
fclose(tipTimesFileID);
fclose(alignFileID);

lastTipTime = max(tipDates)


end

