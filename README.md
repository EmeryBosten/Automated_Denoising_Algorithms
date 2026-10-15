# Automated_Denoising_Algorithms
Matlab codes to "Automated tuning of denoising algorithms in liquid chromatography". 

WS implemtation follows ***"A perfect smoother" by P. Eilers***

SASS implemntation follows ***"Sparsity-Assisted Signal Smoothing (Revisited)" by Ivan Selesnick***

BEADS implementation follows ***"Chromatogram baseline estimation and denoising using sparsity (BEADS)" by Xiaoran Ning, Ivan Selesnick, and Laurent Duval***

demo_script shows a use case on an experimental chromatogram (stored in data) of the different smoother-tuner combinations. 

the experimental chromatogram comes from the paper ***"Baseline correction using adaptive iteratively reweighted penalized least squares" by Zhi-Min Zhang,
Shan Chena and Yi-Zeng Liang*** 

## Demo

```Matlab

%% Demonstration script for "Automated tuning of denoising algorithms for noise removal in liquid chromatography"
% This script gives an example use case of the different smoother-tuner combinations on an example chromatogram 
%% Load data
% data comes from paper 'Baseline correction using adaptive iteratively
% reweighted penalized least squares' by Zhi-Min Zhang, Shan Chena and Yi-Zeng Liang  
% addition of gaussian noise
clear
% add path to data and functions 
addpath('C:\Users\emery\OneDrive\Documents\MATLAB\Automated_Denoising_Algorithms-main\Automated_Denoising_Algorithms-main')

load noise.mat
load chromatograms.mat;

% original chromatogram
x = X(:,5);
% add gaussian noise 
y = x + noise * 1;

% show noisy chromatogram
figure 
plot(y)
hold on 
plot(x)
xlabel('time (sample)','FontSize',15)
ylabel('Absorbance','FontSize',15)
legend('noisy chromatogram','original chromatogram')

%% Initialization of parameter range
% parameter range for Savitzky-Golay filter
SG_ws = 5:2:53;
% parameter range for Whittaker smoother
WS_lambdas = logspace(-3,5,25);
% parameter range for SASS
SASS_lambdas = logspace(-3,5,25);
% parameter range for BEADS
BEADS_alphas = logspace(-3,3,25);
%% Initialization of estimator type
% select estimator type : 'median' or 'mean'
type = 'median';
%% Autocorrelation %%
% autocorrelation function as parameter tuner 
%% Auto-SG
% autocorrelation function in combination with Savitzky-Golay filter
tic
z = autoCorrSG(y,SG_ws,type);
toc
figure 
plot(y,'b')
hold on 
plot(z,'r')
hold on
plot(x)
legend('noisy','estimated ','origin')
title("Auto-SG median estimator")
%% Auto-WS
% autocorrelation function in combination with Whittaker smoother
tic
z = autoCorrWS(y,WS_lambdas,type);
toc
figure 
plot(y,'b')
hold on 
plot(z,'r')
hold on
plot(x)
legend('noisy','estimated','origin')
title("Auto-WS median estimator")
%% Auto-SASS
% autocorrelation function in combination with SASS
tic
z = autoCorrSASS(y,SASS_lambdas,type);
toc
figure 
plot(y,'b')
hold on 
plot(z,'r')
hold on
plot(x)
legend('noisy','estimated','origin')
title("Auto-SASS median estimator")
%% Auto-BEADS
% autocorrelation function in combination with BEADS
tic
[z,b] = autoCorrBEADS(y,BEADS_alphas,type);
toc;
figure 
plot(y,'b')
hold on 
plot(z+b,'r')
hold on
plot(x)
legend('noisy','estimated','origin')
title("Auto-BEADS median estimator")

%% Cross - validation %%
% Cross-validation as parameter tuner 
% choose number of folds
K = 10;
%% CV-SG
% cross-validation in combination with Savitzky-Golay filter
tic
[z]= CVSG(y,SG_ws, K,type);
toc
figure
plot(z)
hold on
plot(y)
hold on
plot(x)
legend('estimated','noisy','origin')
title('CV-SG median estimator')
%% CV-WS
% cross-validation in combination with Whittaker smoother
tic
[z] = CVWS(y,WS_lambdas, K, type);
toc
figure
plot(z)
hold on
plot(y)
hold on
plot(x)
legend('estimated','noisy','origin')
title('CV-WS median estimator')
%% CV-SASS
% cross-validation in combination with SASS
tic
z = CVSASS(y,SASS_lambdas,K,type);
toc
figure
plot(z)
hold on
plot(y)
hold on
plot(x)
legend('estimated','noisy','origin')
title('CV-SASS median estimator')
%% CV-BEADS
% cross-validation in combination with BEADS
tic
[z,b] = CVBEADS(y,BEADS_alphas,K,type);
toc;
figure
plot(z+b)
hold on
plot(y)
hold on
plot(x)
legend('estimated','noisy','origin')
title('CV-BEADS median estimator')
%% V-curve %%
% V-curve as parameter tuner
%% V-Curve-WS
% V-curve in combination with Whittaker smoother
tic
z = VcurveWS(y,WS_lambdas);
toc
figure 
plot(y,'b')
hold on 
plot(z,'r')
hold on
plot(x)
legend('noisy','estimated','origin')
title("V-Curve-WS median estimator")
%% V-Curve-SASS
% V-curve in combination with SASS
tic
z = VcurveSASS(y,SASS_lambdas);
toc
figure 
plot(y,'b')
hold on 
plot(z,'r')
hold on
plot(x)
legend('noisy','estimated','origin')
title("V-Curve-SASS median estimator")
%% LV-Curve-BEADS
% V-curve in combination with BEADS
tic
[z,b] = VcurveBEADS(y,BEADS_alphas);
toc
figure 
plot(y,'b')
hold on 
plot(z+b,'r')
hold on
plot(x)
legend('noisy','estimated','origin')
title('V-Curve-BEADS median estimator')
%% residual Variance %%
% residual variance as parameter tuner 
%% RV-SG
% residual variance in combination with Savitzky-Golay filter
tic
z = resVarSG(y,SG_ws,type);
toc
figure 
plot(y,'b')
hold on 
plot(z,'r')
hold on
plot(x)
legend('noisy','estimated','origin')
title("RV-SG median estimator")
%% RV-WS
% residual variance in combination with Whittaker smoother
tic
z = resVarWS(y,WS_lambdas,type);
toc
figure 
plot(y,'b')
hold on 
plot(z,'r')
hold on
plot(x)
legend('noisy','estimated','origin')
title("RV-WS median estimator")
%% RV-SASS
% residual variance in combination with SASS
tic
z = resVarSASS(y,SASS_lambdas,type);
toc
figure 
plot(y,'b')
hold on 
plot(z,'r')
hold on
plot(x)
legend('noisy','estimated','origin')
title("RV-SASS median estimator")
%% RV-BEADS
% residual variance in combination with BEADS
tic
[z,b] = resVarBEADS(y,BEADS_alphas,type);
toc
figure 
plot(y,'b')
hold on 
plot(z+b,'r')
hold on
plot(x)
legend('noisy','estimated','origin')
title("RV-BEADS median estimator") 
```

## License
see details about the license at "PolyForm NonCommercial 1.0.0.txt"

## Citation

This code is part of the publication "Automated tuning of denoising algorithms for noise removal in chromatograms"

If you use this work, make sure to cite it properly 
Emery Bosten, Peter Van Broeck, Deirdre Cabooter,
Automated tuning of denoising algorithms for noise removal in chromatograms,
Journal of Chromatography A, Volume 1709, 2023, 464360,
https://doi.org/10.1016/j.chroma.2023.464360.
