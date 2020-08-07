function fRho=fisherTransform(ccArray)

fRho = 0.5*(log(1 + ccArray) - log(1.000 - min(ccArray,0.9999)));


