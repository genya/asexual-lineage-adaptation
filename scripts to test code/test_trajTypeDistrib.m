%test trajTypeDistrb()

clear; fclose('all');


%use makeRunFiles.m to generate jobSpecFile
jobSpecFile = 'expDFE_jobSpecs';
timeLimit = 10;
jobInd = 300;

outFile = 'trajTypeDistrib_testOutput';
outFolder = 'trajTypeDistrib_testOutput';
mkdir(outFolder);

trajTypeDistrib(jobInd,timeLimit,outFile,outFolder,jobSpecFile);