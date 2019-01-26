function [ output_args ] = writeBeastXML(filePrefix, xmlTemplate)
% Write simulated data and params to BEAST XML file

if ~exist('filePrefix','var')
  filePrefix = 'fitBD_neutralFitModel_5sites_';
end

% Input template
fid = fopen(xmlTemplate);

% Alignment file
align_file = strcat(filePrefix, '_bdmmAlignment.txt');
align_fid = fopen(align_file);

% Dates file
dates_file = strcat(filePrefix, '_bdmmTipTimes.txt');
dates_fid = fopen(dates_file);

% Types data file
types_file = strcat(filePrefix, '_bdmmTipStates.txt');
types_fid = fopen(types_file);

% Tree data file
tree_file = strcat(filePrefix, '_newick.tre');
tree_fid = fopen(tree_file);

% Output
fout = strcat(filePrefix, '.xml');

% Read all lines into 'xml' cell array
xml = "";
n = 0;
while ~feof(fid)
    n = n+1;
    xml(n) = fgetl(fid);
end
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find start and end of seq alignment
start_idx = find(contains(xml,'<data id="alignment" name="alignment" dataType="binary" spec="beast.evolution.alignment.UnsortedAlignment">')==1);
locs = find(contains(xml,'<sequence id=')==1);
end_idx = locs(end) + 1;

% Splice in new seq alignment
head = xml(1:start_idx);
tail = xml(end_idx:end);
seqs = "";
n = 0;
while ~feof(align_fid)
    n = n+1;
    seqs(n) = fgetl(align_fid);
end
fclose(align_fid);
xml = [head,seqs,tail];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find start and end of date trait block
start_idx = find(contains(xml,'<trait id="dateTrait" spec="beast.evolution.tree.TraitSet" traitname="date-forward">')==1);
locs = find(contains(xml,'<taxa id="TaxonSet" spec="TaxonSet">')==1);
end_idx = locs(end);

% Splice in new date trait black
head = xml(1:start_idx);
tail = xml(end_idx:end);
dates = "";
n = 0;
while ~feof(dates_fid)
    n = n+1;
    dates(n) = fgetl(dates_fid);
end
fclose(dates_fid);
xml = [head,dates,tail];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find start and end of type trait block
start_idx = find(contains(xml,'<trait id="typeTrait" spec="beast.evolution.tree.TraitSet" traitname="type">')==1);
locs = find(contains(xml,'<taxa idref="TaxonSet"/>')==1);
end_idx = locs(end);

% Splice in new date trait black
head = xml(1:start_idx);
tail = xml(end_idx:end);
types = "";
n = 0;
while ~feof(types_fid)
    n = n+1;
    types(n) = fgetl(types_fid);
end
fclose(types_fid);
xml = [head,types,tail];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find start and end of tree block
start_idx = find(contains(xml,'<input name=''newick''>')==1);
locs = find(contains(xml,'</input>')==1);
end_idx = locs(end);

% Splice in new tree
head = xml(1:start_idx);
tail = xml(end_idx:end);
tree = "";
n = 0;
while ~feof(tree_fid)
    n = n+1;
    tree(n) = fgetl(tree_fid);
end
fclose(tree_fid);

% For MultiTypeNexus format
%tree_str = tree{3};
%splits = strsplit(tree_str,'TREE =');
%tree = splits{2};

xml = [head,tree,tail];

% Change log file names
log_idx = find(contains(xml,'<logger fileName="" id="tracelog" logEvery="100" model="@posterior" sanitiseHeaders="true" sort="smart">')==1);
log_file = strcat('<logger fileName="', filePrefix, '.log', '" id="tracelog" logEvery="100" model="@posterior" sanitiseHeaders="true" sort="smart">');
xml(log_idx) = log_file;

% Change stateTrajs file names
%traj_idx = find(contains(xml,'<logger id="trajectLogger" fileName="" logEvery="100">')==1);
%log_file = strcat('<logger id="trajectLogger" fileName="', filePrefix, '_stateTrajs.csv', '" logEvery="100">');
%xml(traj_idx) = log_file;

% Change log file names
%incidence_idx = find(contains(xml,'<logger id="incidenceLogger" fileName="" logEvery="100">')==1);
%log_file = strcat('<logger id="incidenceLogger" fileName="', filePrefix, '_incidence.csv', '" logEvery="100">');
%xml(incidence_idx) = log_file;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write out xml
fid = fopen(fout, 'w');
for n = 1:length(xml)
    fprintf(fid, '%s', xml{n});
    fprintf(fid, '\n');
end
fclose(fid);

end

