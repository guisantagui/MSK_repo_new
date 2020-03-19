basedir = '/Users/santamag/Desktop/GUILLEM/pics/';
cd(basedir)

% get plate file names
swarmtables = dir(fullfile(basedir, 'swarmData_*.xlsx'));
swarmTblName={swarmtables.name};
noTbl = length(swarmTblName);
T = readtable(swarmTblName{1}, 'readVariableNames', true);
for i=2:noTbl
    t = readtable(swarmTblName{i}, 'readVariableNames', true);
    T = [T; t];
end
% fitlm() fitlme() fitglm()

writetable(T, 'allSwarms.csv')