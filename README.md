# Automated_Denoising_Algorithms
Matlab codes to "Automated tuning of denoising algorithms in liquid chromatography". 

WS implemtation follows ***"A perfect smoother" by P. Eilers***

SASS implemntation follows ***"Sparsity-Assisted Signal Smoothing (Revisited)" by Ivan Selesnick***

BEADS implementation follows ***"Chromatogram baseline estimation and denoising using sparsity (BEADS)" by Xiaoran Ning, Ivan Selesnick, and Laurent Duval***

demo_script shows a use case on an experimental chromatogram (stored in data) of the different smoother-tuner combinations. 

the experimental chromatogram comes from the paper ***"Baseline correction using adaptive iteratively reweighted penalized least squares" by Zhi-Min Zhang,
Shan Chena and Yi-Zeng Liang*** 

## Demo


##Load chromatogram
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
```
<img width="1920" height="926" alt="original_chrom" src="https://github.com/user-attachments/assets/40d29abc-edcb-4495-9fc7-8a640d1abea3" />


## Initialize parameters 
```Matlab

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
```

## Run different denoising methods 

```Matlab

%% Autocorrelation %%
% autocorrelation function as parameter tuner 
%% Auto-SG
% autocorrelation function in combination with Savitzky-Golay filter
tic
z = autoCorrSG(y,SG_ws,type);
toc

%% Auto-WS
% autocorrelation function in combination with Whittaker smoother
tic
z = autoCorrWS(y,WS_lambdas,type);
toc

%% Auto-SASS
% autocorrelation function in combination with SASS
tic
z = autoCorrSASS(y,SASS_lambdas,type);
toc

%% Auto-BEADS
% autocorrelation function in combination with BEADS
tic
[z,b] = autoCorrBEADS(y,BEADS_alphas,type);
toc;

%% Cross - validation %%
% Cross-validation as parameter tuner 
% choose number of folds
K = 10;
%% CV-SG
% cross-validation in combination with Savitzky-Golay filter
tic
[z]= CVSG(y,SG_ws, K,type);
toc

%% CV-WS
% cross-validation in combination with Whittaker smoother
tic
[z] = CVWS(y,WS_lambdas, K, type);
toc

%% CV-SASS
% cross-validation in combination with SASS
tic
z = CVSASS(y,SASS_lambdas,K,type);
toc

%% CV-BEADS
% cross-validation in combination with BEADS
tic
[z,b] = CVBEADS(y,BEADS_alphas,K,type);
toc;

%% V-curve %%
% V-curve as parameter tuner
%% V-Curve-WS
% V-curve in combination with Whittaker smoother
tic
z = VcurveWS(y,WS_lambdas);
toc

%% V-Curve-SASS
% V-curve in combination with SASS
tic
z = VcurveSASS(y,SASS_lambdas);
toc

%% LV-Curve-BEADS
% V-curve in combination with BEADS
tic
[z,b] = VcurveBEADS(y,BEADS_alphas);
toc

%% residual Variance %%
% residual variance as parameter tuner 
%% RV-SG
% residual variance in combination with Savitzky-Golay filter
tic
z = resVarSG(y,SG_ws,type);
toc

%% RV-WS
% residual variance in combination with Whittaker smoother
tic
z = resVarWS(y,WS_lambdas,type);
toc

%% RV-SASS
% residual variance in combination with SASS
tic
z = resVarSASS(y,SASS_lambdas,type);
toc

%% RV-BEADS
% residual variance in combination with BEADS
tic
[z,b] = resVarBEADS(y,BEADS_alphas,type);
toc
<img width="1920" height="926" alt="original_chrom" src="https://github.com/user-attachments/assets/b859c10b-8f28-4bb0-a03e-2f3ad22f0caf" />

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
