function [] = testcontourstereogram()
 commandwindow
keyboard
%% Simple test
close all; clear all; clc
eulerAngles = [35,16,280; 173,163,023; 267,016,163; 000,173,123]*pi/180;
mineral = 'olivine';
contourpolefigures(eulerAngles,mineral,'Gaussian',1,'lower',51)
%% Test Olivine
close all; clear all; clc
% Read in euler angle file
fileName = 'drex_simpleShear.txt';
fid = fopen(fileName);
angles = textscan(fid,'%f %f %f %*[^\r\n]','HeaderLines',1);
eulerAngles = [angles{1},angles{2},angles{3}]*pi/180;
clear angles;

mineral = 'olivine';
method = 'Gaussian';
ColorOpts.colorRamp0 = coolwarm(256);
ColorOpts.nTics      = 3;
ColorOpts.centerVal  = 1;
ColorOpts.factor     = 0.5;

hFig = contourpolefigures(eulerAngles,mineral,method,[],[],[],ColorOpts);


%% Test Quartz
clc
close all  
% fileName = 'test-quartz.txt';
% fileName = 'MT-07-76.txt';
% fileName = 'MT1a-out.txt';
fileName = 'MT6-out.txt';
fid = fopen(fileName);
angles = textscan(fid,'%f %f %f %*[^\r\n]','HeaderLines',1);
eulerAngles = [angles{1},angles{2},angles{3}]*pi/180;
clear angles


mineral = 'quartz';
method = 'Gaussian';
ColorOpts.colorRamp0 = coolwarm(256);
ColorOpts.nTics      = 3;
ColorOpts.centerVal  = 1;
ColorOpts.factor     = 0.5;

hFig = contourpolefigures(eulerAngles,mineral,method,[],[],151,ColorOpts);


end

