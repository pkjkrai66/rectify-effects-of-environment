% This script evaluates the performance of multiple functions by regressing
% and calculate the RMSE and displays in table format
clear;
addpath('H:\Academics\Matlab Codes\Matlab Support Functions')

%% Only this section need to be changed                                 |
% load('DWL_T_Internal_Temp_Compensation.mat')                  %       |
% DataSet = DWLTData;                                           %       |
                                                                %       |
load('UWL_NTC_CombinedData.mat')                                %       |
DataSet = UWL_NTC_DataCombined;                                 %       |
    % Remove measurements where DataSet == NaN                          |
notargetIndices = isnan(DataSet.Target);                        %       |
DataSet         = DataSet(~notargetIndices, :);                 %       |
timeVector      = DataSet.Timestamp;                            %       |
xVector         = DataSet.mode;   % Measurement from instrument         |
temperatureVec  = DataSet.Temp;     % Ambient Temperature               |
TargetVec       = DataSet.Target;   % Target = true measurement         |
                                                                %       |
X = [xVector(:), temperatureVec(:)];                            %       |
Y = TargetVec;                                                  %       |
                                                                %       |
% ______________________________________________________________________|

clear notargetIndices 
%% Declaring Functions ------------------
funlist = {'b(1) + b(2).*X(:,1) + b(3)*log(X(:,2)).*X(:,1)';...
        'b(1) + b(2).*X(:,1) + sign(X(:,2)-23).*b(3).*log(X(:,2)).*X(:,1)'};
        % list of all the functions that need to be evaluated
nfun    = size(funlist,1);     % No of functions
hndl1   = repmat({'@(b, X)'}, [nfun,1]); 
strfun  = strcat(hndl1, funlist);
for i=1:nfun
   fun{i,1} = str2func(strfun{i});
end
% --------------------------------------------------
% Or rather a simple way ----------------
% fun(:,1) = {@(b, X)(b(1) + b(2).*X(:,1) + b(3)*log(X(:,2)).*X(:,1));...
%     @(b, X)(b(1) + b(2).*X(:,1) + sign(X(:,2)-23).*b(3).*log(X(:,2)).*X(:,1))};
clear funlist hndl1 strfun 
beta0(:,1) = {[0.03, 1, 10];...
            [0.03, 1, 0.10]};
%% Regression
for i = 1:nfun
    betaE{i} =      nlinfit(X, Y, fun{i,1}, beta0{i});
    Yestimate{i} =  feval(fun{i,1}, betaE{i}, X);
end
clr = jet(2+nfun);
mark    = {'o','sq','x','.','+','pentagram'};
markSz  = [  2,   3,  3,  6,  3,         2];
lineWd  = [0.5,   1,  1,  2,  1,         1];

clear hndl1 betaE X Y
%% Variables
Targets = sort(unique(TargetVec))'; % Targets= distance for which measurements are taken

%% Calculating R2 and RMSE

% Calculating distance using Manufacturer Recomended Equation
% Follow the link for details- 
%   https://www.maxbotix.com/documents/Temperature_Compensation.pdf
TOF = xVector*(147*10^(-6)/25.4);   % Convert distance to time Of Flight
Tc = temperatureVec;                % Ambient Temperature in Centigrade
Dm = TOF.*20.05.*sqrt(Tc+273.15)*1000/2;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

RMSE = RootMeanSqrErr(TargetVec, xVector);
i = 0;
% Calculating RMSE for a given target `TargetC`
for targetValue = Targets
    I = (targetValue == TargetVec);   % Indices >> I 
    i = i+1;
    
    RMSE_obser(i,1) = RootMeanSqrErr(TargetVec(I), xVector(I));
    RMSE_Manufacturer(i,1) = RootMeanSqrErr(TargetVec(I), Dm(I));
    targetStr{i} = sprintf('t=%d mm', targetValue); 
    for j = 1:nfun
        RMSE_fun(i,j) = RootMeanSqrErr(TargetVec(I), Yestimate{j}(I));
        funStr{j} = sprintf('fun%d', j);
    end
    
end
tb0 = array2table(RMSE_fun, 'VariableNames', funStr);
tb1 = table(RMSE_obser, RMSE_Manufacturer, tb0, 'VariableNames',...
   {'RMSE_observations', 'RMSE_Manufacturer', 'RMSE_modelled'}, 'RowNames', targetStr);
disp(tb1)

%% Plotting 
% -----------------------

% Figure
f1 = figure;
axf1 = axes();
hold on; box on;
i = 0;
for targetValue = Targets
    I = (targetValue == TargetVec);
    i = i+1;
    p1{i} = plot(timeVector(I), TargetVec(I), 'k', 'linewidth', 1);  % Target
    p2{i} = plot(timeVector(I), xVector(I), 'Color', clr(1, :),...
            'Marker', mark{1}, 'MarkerSize', markSz(1), ...
            'lineStyle', 'none', 'LineWidth', lineWd(1));
    p3{i} = plot(timeVector(I), Dm(I), 'Color', clr(2, :),...
            'Marker', mark{2}, 'MarkerSize', markSz(2), ...
            'lineStyle', 'none', 'LineWidth', lineWd(2));
    % Plotting for the fitted functions
    for j = 1:nfun
        hndlfun1(i,j) = plot(timeVector(I), Yestimate{j}(I), 'Color', clr(2+j,:),...
            'Marker', mark{2+j}, 'MarkerSize', markSz(2+j), ...
            'lineStyle', 'none', 'LineWidth', lineWd(2+j)); % Estimate
        temp_st = char(fun{j});
        tempCell{j} = temp_st(7:end);
    end
end
legend([p1{i}, p2{end}, p3{end}, hndlfun1(i,:)], ...
      {'Target', 'Instrument reading', 'Manufacturer Cal Eqn.', tempCell{:}});
xlabel('time'); ylabel('distance (mm)')

% Plot corrpsponding to a target value
f2 = figure;
i = 0;
for targetValue = Targets
    I = (targetValue == TargetVec);
    i = i + 1;
    axf2(i) = subplot(2, 3, i); 
    hold on; box on; 
    ps1 = plot(timeVector(I), TargetVec(I), 'k', 'linewidth',1);  % Target
    ps2 = copyobj(p2{i}, axf2(i));
    ps3 = copyobj(p3{i}, axf2(i));
    for j = 1:nfun
        hndlfun(j) = copyobj(hndlfun1(i,j), axf2(i)); % Estimate
        temp_st = sprintf('fun%d Est (RMSE %.0f)', [j, RMSE_fun(i,j)]);
        tempCell{j} = temp_st; 
    end
    xlabel('time'); ylabel('distance (mm)');
    title(sprintf('Target = %d', targetValue));
end
legend([ps1, ps2, ps3, hndlfun(:)'], {'Target', ...
        sprintf('Instrument reading (RMSE: %0.0f)', RMSE_obser(i)), ...
        sprintf('Manufactur Cal Eqn.(RMSE: %0.0f)', RMSE_Manufacturer(i)), ...
        tempCell{:}})
clear p1 p2 axf2 ps1 ps2 
%% plotting observation(y) v/s temperature(x) for a given target distance
f3 = figure;
i = 0;
for targetValue = Targets
    I = (targetValue == TargetVec);
    i = i+1;
    axf3(i) = subplot(2, 3, i); 
    hold on; box on; 
    title(['target = ', num2str(targetValue), ' mm'])
    p21 = plot(temperatureVec(I), TargetVec(I), 'k', 'linewidth',1);
    p22 = plot(temperatureVec(I), xVector(I), 'Color', clr(1, :),...
            'Marker', mark{1}, 'MarkerSize', markSz(1), ...
            'lineStyle', 'none', 'LineWidth', lineWd(1));
    p23 = plot(temperatureVec(I), Dm(I), 'Color', clr(2, :),...
            'Marker', mark{2}, 'MarkerSize', markSz(2), ...
            'lineStyle', 'none', 'LineWidth', lineWd(2));
    for j = 1:nfun
        hndlfun2(j) = plot(temperatureVec(I), Yestimate{j}(I), 'Color', clr(2+j, :),...
            'Marker', mark{2+j}, 'MarkerSize', markSz(2+j), ...
            'lineStyle', 'none', 'LineWidth', lineWd(2+j)); % Estimate
    end
end
legend([p21, p22, p23, hndlfun2(:)'], {'Target', ...
        sprintf('Instrument reading (RMSE: %0.0f)', RMSE_obser(i)), ...
        sprintf('Manufactur Cal Eqn.(RMSE: %0.0f)', RMSE_Manufacturer(i)), ...
        tempCell{:}})


