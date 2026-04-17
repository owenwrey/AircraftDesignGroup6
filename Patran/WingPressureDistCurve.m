clear
close all
clc

n = 50;
m = 10;
b = 2;
y = linspace(0,1,20);
x = linspace(0.1,0.85,n);


spanwisePressureDistributionNormalized = sqrt(1 - (2*y/b).^2);

upperSurfacePressureDistData = readtable("UpperSurfacePressureDist.csv");
lowerSurfacePressureDistData = readtable("LowerSurfacePressureDist.csv");

chordwisePressureDistUpperSurfaceNorm = interp1(upperSurfacePressureDistData.x,upperSurfacePressureDistData.y,x,"makima","extrap");
chordwisePressureDistLowerSurfaceNorm = interp1(lowerSurfacePressureDistData.x,lowerSurfacePressureDistData.y,x,"makima","extrap");

PatranCoords = table( ...
[5.3465301e-06; 5.9864346e-06; 8.1874568e-06; 7.4365398e-06; 1.1092272e-05; ...
 1.0311406e-05; 1.4049880e-05; 1.3146931e-05; 1.7043434e-05; 1.5962978e-05; ...
 1.7369792e-05], ...

[9.6873732; -4.9532495; 10.693595; -5.1280718; 10.23794; ...
 -4.8978791; 8.5760593; -3.7667503; 6.0958991; -2.1919506; ...
 -1.3612754], ...

[-65.999138; -65.999802; -32.999626; -49.499943; -0.00024221298; ...
 -16.50021; 32.999905; 16.499973; 66.002205; 49.499207; ...
 66.006599], ...

'VariableNames', {'x','y','z'} );
xReal = linspace(max(PatranCoords.z),min(PatranCoords.z),30);
yReal = linspace(max(PatranCoords.x),min(PatranCoords.x),20);
LECoords = linspace(PatranCoords.z(1), PatranCoords.z(2),20);
TECoords = linspace(PatranCoords.z(3), PatranCoords.z(4),20);


for i=1:20
ChordCoords = linspace(LECoords(i), TECoords(i),n);

UpperPressureDist(i,:) = interp1(ChordCoords,chordwisePressureDistUpperSurfaceNorm,xReal,"nearest","extrap");
LowerPressureDist(i,:) = interp1(ChordCoords,chordwisePressureDistLowerSurfaceNorm,xReal,"nearest","extrap");
end

UpperPressureDist = UpperPressureDist.*repmat(spanwisePressureDistributionNormalized',1,30);
LowerPressureDist = LowerPressureDist.*repmat(spanwisePressureDistributionNormalized',1,30);

l = 1;
for m=1:length(xReal)
    for n = 1:length(yReal)
        UpperPressureField(l,:) = [yReal(n),xReal(m),UpperPressureDist(n,m)];
        LowerPressureField(l,:) = [yReal(n),xReal(m),LowerPressureDist(n,m)];
        l = l+1;
    end
end


writematrix(UpperPressureField,"UpperPressureField.csv","Delimiter",'comma')
writematrix(LowerPressureField,"LowerPressureField.csv","Delimiter",'comma')