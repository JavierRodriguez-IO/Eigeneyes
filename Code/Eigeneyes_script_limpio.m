% Eigeneyes: A statistical model of whole-eye geometry and its applications
% Copyright (C) 2026 Estela Vadillo, Álvaro de la Peña,
% Javier Rodríguez-Sánchez, Alberto de Castro, Susana Marcos,
% Eduardo Martínez-Enríquez
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

%% EIGENEYES - MODELO FINAL Y VALIDACION CRUZADA 10 x 10
clear; close all; clc;
rng(42);

%% CONFIGURACION
nFolds = 10;
nRepeats = 10;
Ks = 1:10;
K = numel(Ks);
nK = numel(Ks);

%% DATOS
scriptFolder = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptFolder,'functions'));
resultsFolder = fullfile(scriptFolder,'..','Results');
if ~isfolder(resultsFolder), mkdir(resultsFolder); end

dataFolder = fullfile(scriptFolder,'..','Data');
dataFile = fullfile(dataFolder,'eigeneyes_complete_dataset_all_right.mat');
if ~isfile(dataFile)
    error('No se encuentra el dataset: %s\nColoque eigeneyes_complete_dataset_all_right.mat en la carpeta Data.', dataFile);
end

load(dataFile);

% Se esperan dataZ_eigeneyes, dataX_eigeneyes, dataY_eigeneyes, Age e ID.
dataZ_eigeneyes_pre = dataZ_eigeneyes;
dataX_eigeneyes_pre = dataX_eigeneyes;
dataY_eigeneyes_pre = dataY_eigeneyes;
age_vector = round(Age(:));
Num_eyes = size(dataZ_eigeneyes_pre,2);
ID = string(ID(:));
assert(numel(ID)==Num_eyes,'Debe existir un ID de sujeto por ojo.');

% Aplicar solo si los datos no estan ya registrados en retina.
dataZ_eigeneyes_pre = eigeneyes_image_registration('retina', dataZ_eigeneyes_pre, Num_eyes);

%% EIGENEYES DEFINITIVOS (TODOS LOS OJOS)
Kfinal = max(Ks);
[RMS,kp5,coefs_correlation,polynomials_third_order, ...
 Reconstructed_eyes,signals,tv,mn,PC,cum_explained] = ...
    eigeneyes_SVD_cartesian( ...
    dataZ_eigeneyes_pre,dataX_eigeneyes_pre,dataY_eigeneyes_pre, ...
    0.99,Num_eyes,Kfinal,age_vector);

save(fullfile(resultsFolder,'Eigeneyes_modelo_final.mat'),'mn','PC','signals','cum_explained', ...
    'RMS','Reconstructed_eyes','kp5','ID','-v7.3');

%% VALIDACION CRUZADA AGRUPADA POR SUJETO
X = double(dataZ_eigeneyes_pre.');  % ojos x puntos
[uniqueSubjects,~,subjectIdx] = unique(ID);
nSubjects = numel(uniqueSubjects);
assert(nSubjects>=nFolds,'No hay suficientes sujetos para 10 folds.');
assert(all(isfinite(X(:))),'Hay valores NaN o Inf.');

RMSE_CV = nan(nRepeats,nFolds,nK);
PVAR_CV = nan(nRepeats,nFolds,nK);
RMSE_eye_CV = nan(Num_eyes,nK,nRepeats);
subjectFold_CV = zeros(nSubjects,nRepeats);

for rep = 1:nRepeats
    cv = cvpartition(nSubjects,'KFold',nFolds);
    fprintf('Repetition %d/%d\n',rep,nRepeats);

    for fold = 1:nFolds
        trainSubjects = training(cv,fold);
        testSubjects = test(cv,fold);
        trainIdx = trainSubjects(subjectIdx);
        testIdx = testSubjects(subjectIdx);
        subjectFold_CV(testSubjects,rep) = fold;

        Xtrain = X(trainIdx,:);
        Xtest = X(testIdx,:);
        mu = mean(Xtrain,1);
        Xtrain0 = Xtrain-mu;
        Xtest0 = Xtest-mu;

        % PCA del fold: solo entrenamiento, sin fuga de informacion.
        [~,S,V] = svd(Xtrain0,'econ');
        eigval = diag(S).^2;
        pvar = 100*cumsum(eigval)/sum(eigval);
        assert(max(Ks)<=size(V,2),'K supera el rango disponible.');

        for j = 1:nK
            K = Ks(j);
            Vk = V(:,1:K);
            Xrec = (Xtest0*Vk)*Vk' + mu;
            rmseEye = sqrt(mean((Xtest-Xrec).^2,2));

            RMSE_CV(rep, fold, j) = mean(rmseEye);
            PVAR_CV(rep,fold,j) = pvar(K);
            RMSE_eye_CV(testIdx,j,rep) = rmseEye;
        end
    end
end

%% RESULTADOS
RMSE_all = reshape(RMSE_CV,[],nK);
PVAR_all = reshape(PVAR_CV,[],nK);
RMS_vector = mean(RMSE_all,1);
RMS_std_vector = std(RMSE_all,0,1);
PV_vector = mean(PVAR_all,1);
PV_std_vector = std(PVAR_all,0,1);

Results = table(Ks(:),RMS_vector(:),RMS_std_vector(:), ...
    PV_vector(:),PV_std_vector(:), ...
    'VariableNames',{'K','RMSE_mean_mm','RMSE_std_mm', ...
    'PVAR_mean_percent','PVAR_std_percent'});
disp(Results);

save(fullfile(resultsFolder,'Eigeneyes_10x10_CV.mat'),'Ks','nFolds','nRepeats', ...
    'RMSE_CV','PVAR_CV','RMSE_eye_CV','subjectFold_CV', ...
    'RMS_vector','RMS_std_vector','PV_vector','PV_std_vector', ...
    'Results','-v7.3');
writetable(Results,fullfile(resultsFolder,'Eigeneyes_10x10_CV_results.csv'));

%% GRAFICA
figure('Color','w','Position',[100 100 1200 750]);
yyaxis left
h1 = errorbar(Ks,RMS_vector,RMS_std_vector,'o-','LineWidth',2);
ylabel('RMSE (mm)');
yyaxis right
h2 = errorbar(Ks,PV_vector,PV_std_vector,'s--','LineWidth',2);
ylabel('PVar (%)');
xlabel('Number of Eigeneyes (K)');
xticks(Ks);
legend([h1 h2],{'RMSE','PVar'},'Location','best');
set(gca,'FontSize',18,'LineWidth',1.5);
box on;
exportgraphics(gcf,fullfile(resultsFolder,'Eigeneyes_10x10_CV_RMSE_PVAR.png'),'Resolution',300);
