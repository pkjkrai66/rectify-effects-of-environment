% The purpose of this script is to test the performance of 
%   logistics function on data
clear;
load('DWL_T_Internal_Temp_Compensation.mat')
addpath('H:\Academics\Matlab Codes\Matlab Support Functions')

DataSet = DWLTData;

timeVector  = DataSet.Timestamp;
xVector     = DataSet.mode;         % Measurement from instrument
temperatureVec  = DataSet.TempExt;  % Ambient Temperature 
TargetVec   = DataSet.TargetCorrect;% Target = true measurement 

x = 0:0.1:40;
L = 50;
x0 = 23;
k = 0.55;

f4 = figure;

%% Variables
Targets = sort(unique(TargetVec))'; % Targets= distance for which measurements are taken

i = 0;
for targetValue = Targets
    I = (targetValue == TargetVec);
    i = i+1;
    axf4(i) = subplot(2, 3, i); 
    f1zx1 = plot(temperatureVec(I), xVector(I), 'r.');
    hold on;

    fun_logistic = targetValue + (L./(1+exp(-k*(x-x0))) - L/2);

    f1zx2 = plot(x, fun_logistic, 'k');
    myHline(targetValue); myVline(x0);
    xlabel('temperature')
    ylabel('distance (mm)')
    title(['target = ', num2str(targetValue)])

end

