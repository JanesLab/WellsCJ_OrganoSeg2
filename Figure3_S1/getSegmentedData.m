load('goldStandardComp_Original.mat');
signatures.borten = signatureStruct;
ccs.borten = hiddenProps.CC;

load('goldStandardComp_EdgeCorrect_Watershed.mat');
signatures.edge = signatureStruct;
ccs.edge = hiddenProps.CC;

