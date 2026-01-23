function T2W = TakeoffConstraintAnalysis(TOP,sigma,Cl,W2S)
T2W = W2S.*0.02088./TOP./sigma./Cl;
end