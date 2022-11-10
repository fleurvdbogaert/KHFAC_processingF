%% General Information
% Author: Fleur van den Bogaert, fleurvandenbogaert@live.nl
% Technical Medicine Internship September 2022 - November 2022 

% Previous authors: 
    % Max Ligtenberg, max.ligtenberg@outlook.com, June - Aug 2020
    % Bart Formsma, bartformsma@hotmail.com, Sep - Nov 2020
    % Anne Meester, annemeester95@hotmail.com, Sep - Nov 2021
    % Sabine Josemans, shjosemans@gmail.com, Nov 2021 - Feb 2022
    % Yasmin Ben Azouz, yasminbenazouz@hotmail.com, May 2022 - August 2022 
    
% ErasmusMC, dept. Urology, group: functional Urology lab
% Edited in MATLAB R2021b

close all; clear; clc;
%% Variables
FS_S = 60000;        % [Hz] sampling frequency stimulation
FS_P = 60000;        % [Hz] sampling frequency pressure 
%% Step 1. 
[STIM,PRES,FILE,DATA] = loadModify(FS_P,FS_S) ;
%% Step 2a. 
[INT_STIM] = stimDetection(STIM) ; 
%% Step 2b. 
[INT_PRES] = contDetection(PRES) ; 
%% Step 3a. 
[CHCK_STIM] = manualCheck(INT_STIM,'stimulation') ; 
%% Step 3b. 
[CHCK_PRES] = manualCheck(INT_PRES,'pressure') ;
%% Step 4. 
[T,Tav_pres] = calcExport(CHCK_STIM,CHCK_PRES,FILE,DATA,FS_S) ; 