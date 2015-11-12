% The directory and file names of amplitude, real and imaginary data
clear all; clc;
global FPath1
RootPath = 'C:\Research\Data\GaAs 4QW (Anharmonic oscillator)\2015_04_08\';
FPath1 = strcat(RootPath,'2D6\Output\');    % Co-circular

% gAbsFreqData = 'gAbsFreqChopped.dat';
% gEmiFreqData = 'gEmiFreqChopped.dat';

global NPtsExD1 NPtsExX1 ExPkX1 ExPkY1...
    Size ExFlagD1 ExFlagX1 ToFit meV2Hz
global PlotSimSpec OutPath

meV2Hz = 241.79895E9*2*pi;      % hbar
ChkEgyRange = 1;    % Also use to plot and save measured spectra
ChkCutSpec = 0;
ShowSlices = 0;
PlotSimSpec = 0;

Power = 153;    % Per beam in uW
if ~isdir(strcat(RootPath,'AO_fits_CoCir_FixDipole\Pow_',num2str(Power),'uW\'))
    mkdir(RootPath,strcat('AO_fits_CoCir_FixDipole\Pow_',num2str(Power),'uW'));
end
OutPath = strcat(RootPath,'AO_fits_CoCir_FixDipole\Pow_',num2str(Power),'uW\');

%% Read frequency mesh matrices
gAbsFreqData = 'gAbsFreq.dat';
gEmiFreqData = 'gEmiFreq.dat';
EmiFMesh1 = dlmread(strcat(FPath1,gEmiFreqData), '\t');
AbsFMesh1 = dlmread(strcat(FPath1,gAbsFreqData), '\t');
MDim1 = [EmiFMesh1(1) EmiFMesh1(end) -EmiFMesh1(1) -EmiFMesh1(end)];
MSize1 = size(EmiFMesh1);
EmiFStep1 = (EmiFMesh1(end) - EmiFMesh1(1))/(MSize1(2)-1);
AbsFStep1 = (AbsFMesh1(end) - AbsFMesh1(1))/(MSize1(1)-1);
EmiAxis1 = zeros(1,MSize1(2));
AbsAxis1 = zeros(1,MSize1(1));
for j = 1 : MSize1(2)
    EmiAxis1(j) = EmiFMesh1(1) + EmiFStep1*(j-1);
end
for j = 1 : MSize1(1)
    AbsAxis1(j) = AbsFMesh1(1) + AbsFStep1*(j-1);
end

EmiFMesh1 = repmat(EmiAxis1,MSize1(1),1);
AbsFMesh1 = repmat(AbsAxis1',1,MSize1(2));

%% Read the 2D spectra
AbsFName1 = 'MAmpl1.dat';
RealFName1 = 'MReal1.dat';
AbsMData1 = dlmread(strcat(FPath1,AbsFName1), '\t');
RealMData1 = dlmread(strcat(FPath1,RealFName1), '\t');
VMax1 = max(max(AbsMData1));

%% Select the relevant region
FMax = 1548;
FMin = 1540;

[Ra1,~] = find(AbsFMesh1 < -FMax, 1, 'last');
[Rb1,~] = find(AbsFMesh1 > -FMin, 1, 'first');
[~,Ca1] = find(EmiFMesh1 < FMin, 1, 'last');
[~,Cb1] = find(EmiFMesh1 > FMax, 1, 'first');

EmiFMeshCut1 = EmiFMesh1(Ra1:Rb1,Ca1:Cb1);
AbsFMeshCut1 = AbsFMesh1(Ra1:Rb1,Ca1:Cb1);
AbsMDataCut1 = AbsMData1(Ra1:Rb1,Ca1:Cb1);
RealMDataCut1 = RealMData1(Ra1:Rb1,Ca1:Cb1);

if ChkEgyRange
    Fig1 = figure(1);
    set(Fig1,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    subplot(121);
    contourf(EmiFMeshCut1, AbsFMeshCut1, AbsMDataCut1, linspace(0,VMax1,50),...
        'LineStyle','none');
    title('Co-cir Spectrum:Abs')
    axis square;
    line([EmiFMeshCut1(1) EmiFMeshCut1(end)], [AbsFMeshCut1(end) AbsFMeshCut1(1)],...
       'LineStyle', '--', 'Color', [1 1 1],'LineWidth',1.5);
    subplot(122);
    contourf(EmiFMeshCut1, AbsFMeshCut1, RealMDataCut1, linspace(-VMax1,VMax1,50),...
        'LineStyle','none');
    title('Co-cir Spectrum:Real')
    axis square;
    line([EmiFMeshCut1(1) EmiFMeshCut1(end)], [AbsFMeshCut1(end) AbsFMeshCut1(1)],...
       'LineStyle', '--', 'Color', [1 1 1],'LineWidth',1.5);
    saveas(gcf,strcat(OutPath,'MeasuredSpectra'),'emf');
    
    ContButton = questdlg('Is the energy range correct?');
    if strcmp(ContButton,'Yes')
    else
        return;
    end
end

% Resize 2D spectra axes and spectra
Size = 512;
% EmiFStep = (EmiFMeshCut1(end) - EmiFMeshCut1(1))/(Size-1);
% AbsFStep = (AbsFMeshCut1(end) - AbsFMeshCut1(1))/(Size-1);
EmiFStep = (FMax - FMin)/(Size-1);
AbsFStep = EmiFStep;
EmiAxis = zeros(1,Size);
AbsAxis = zeros(1,Size);
for j = 1 : Size
    EmiAxis(j) = FMin + EmiFStep*(j-1);
    AbsAxis(j) = -FMax + AbsFStep*(j-1);
end

global EmiFMesh AbsFMesh
EmiFMesh = repmat(EmiAxis,Size,1);
AbsFMesh = repmat(AbsAxis',1,Size);
% EmiFMesh = zeros(Size,Size);
% AbsFMesh = zeros(Size,Size);
% for j = 1:Size
%     EmiFMesh(j,:) = EmiAxis;
%     AbsFMesh(:,j) = AbsAxis';
% end

AbsMData1 = interp2(EmiFMeshCut1, AbsFMeshCut1, AbsMDataCut1, EmiFMesh,...
    AbsFMesh);
RealMData1 = interp2(EmiFMeshCut1, AbsFMeshCut1, RealMDataCut1, EmiFMesh,...
    AbsFMesh);

% Save measured 2D spectrum data and frequency meshes
dlmwrite(strcat(OutPath,'EmiFMesh.dat'),EmiFMesh,'Delimiter','\t','Precision',6);
dlmwrite(strcat(OutPath,'AbsFMesh.dat'),AbsFMesh,'Delimiter','\t','Precision',6);
dlmwrite(strcat(OutPath,'CoCirAbsSpecData.dat'),AbsMData1,'Delimiter','\t','Precision',6);
dlmwrite(strcat(OutPath,'CoCirRealSpecData.dat'),RealMData1,'Delimiter','\t','Precision',6);

if ChkCutSpec
    figure(2);
    subplot(121);
    contourf(EmiFMesh, AbsFMesh, AbsMData1, linspace(0,VMax1,50),...
        'LineStyle','none');
    axis square;
%     line([EmiFMesh(1) EmiFMesh(end)], [AbsFMesh(end) AbsFMesh(1)]);
%     line([EmiFMesh(1) EmiFMesh(end)], [AbsFMesh(1) AbsFMesh(end)]);
    
    subplot(122);
    contourf(EmiFMesh, AbsFMesh, RealMData1, linspace(-VMax1,VMax1,50),...
        'LineStyle','none');
    axis square;
%     line([EmiFMesh(1) EmiFMesh(end)], [AbsFMesh(end) AbsFMesh(1)]);
%     line([EmiFMesh(1) EmiFMesh(end)], [AbsFMesh(1) AbsFMesh(end)]);
    
    ContButton = questdlg('Continue?');
    if strcmp(ContButton,'Yes')
    else
        return;
    end
end

%% Define mesh grids in t and tau
global tMesh tauMesh
MSize = size(EmiFMesh);
tAxis = zeros(1,MSize(2));
tauAxis = zeros(MSize(1),1);
tMesh = zeros(MSize);
tauMesh = tMesh;
tStep = 2*pi/(meV2Hz*(EmiFMesh(1,MSize(2))-EmiFMesh(1,1)));
tauStep = 2*pi/(meV2Hz*(AbsFMesh(MSize(1),1)-AbsFMesh(1,1)));
for j = 1:MSize(2);
    tAxis(1,j) = tStep*(j-1);
end
for j = 1:MSize(1);
    tauAxis(j,1) = tauStep*(j-1);
    tMesh(j,:) = tAxis(:);
end
for j = 1:MSize(2);
    tauMesh(:,j) = tauAxis(:);
end

%% Take relevant slices from the spectra
% Find position of max in absolute spectrum for co-circular polarization
VMax1 = max(max(AbsMData1));
[ExPkY1, ExPkX1] = find(AbsMData1==VMax1);

% Horizontal slices
%Co-circular
AbsHSlice1 = AbsMData1(ExPkY1,:);
RealHSlice1 = RealMData1(ExPkY1,:);

% Vertical slices
% Co-circular
AbsVSlice1 = AbsMData1(:,ExPkX1);
RealVSlice1 = RealMData1(:,ExPkX1);

if ShowSlices
    figure(3);
    subplot(121);
    contourf(EmiFMesh, AbsFMesh, AbsMData1, linspace(0,VMax1,50),...
        'LineStyle','none');
    axis square;
    line([EmiFMesh(1) EmiFMesh(end)], [AbsFMesh(ExPkY1,1) AbsFMesh(ExPkY1,1)]);
    line([EmiFMesh(1,ExPkX1) EmiFMesh(1,ExPkX1)], [AbsFMesh(1) AbsFMesh(end)]);
    
    subplot(122)
    contourf(EmiFMesh, AbsFMesh, RealMData1, linspace(-VMax1,VMax1,50),...
        'LineStyle','none');
    axis square;
    line([EmiFMesh(1) EmiFMesh(end)], [AbsFMesh(ExPkY1,1) AbsFMesh(ExPkY1,1)]);
    line([EmiFMesh(1,ExPkX1) EmiFMesh(1,ExPkX1)], [AbsFMesh(1) AbsFMesh(end)]);
    
    ContButton = questdlg('Continue?');
    if strcmp(ContButton,'Yes')
    else
        return;
    end
end

% Alignment of spectra from simulation and experiment
% Co-circular
PkAbsF1 = AbsFMesh(ExPkY1,ExPkX1);
PkEmiF1 = EmiFMesh(ExPkY1,ExPkX1);
YCor = (PkAbsF1 - (AbsFMesh(1)+AbsFMesh(end))/2);
XCor = (PkEmiF1 - (EmiFMesh(1)+EmiFMesh(end))/2);
XCor1 = XCor*meV2Hz;
YCor1 = YCor*meV2Hz;

% Diagonal slices
% Co-circular
NPtsExD1 = ExPkX1 + ExPkY1 -1;
if NPtsExD1 <= Size
    ExDiagSlice1 = zeros(1,NPtsExD1);
    for j = 1:NPtsExD1
        ExDiagSlice1(j) = AbsMData1(NPtsExD1+1-j,j);
    end
    ExDiagF1 = EmiAxis(1:NPtsExD1);
    ExFlagD1 = 1;
else
    NPtsExD1 = 2*Size - NPtsExD1;
    ExDiagSlice1 = zeros(1,NPtsExD1);
    for j = 1:NPtsExD1
        ExDiagSlice1(j) = AbsMData1(Size+1-j,Size-NPtsExD1+j);
    end
    ExDiagF1 = EmiAxis(Size - NPtsExD1 +1:Size);
    ExFlagD1 = 0;
end

% X-Diagonal slices
% Co-circular
ExDiff1 = ExPkY1 - ExPkX1;
NPtsExX1 = Size - abs(ExDiff1);
ExXDiagSlice1 = zeros(1,NPtsExX1);
ExRealXDiagSlice1 = zeros(1,NPtsExX1);
if ExDiff1 <= 0
    for j = 1 : NPtsExX1
        ExXDiagSlice1(j) = AbsMData1(j,Size-NPtsExX1+j);
        ExRealXDiagSlice1(j) = RealMData1(j,Size-NPtsExX1+j);
    end
    ExXDiagF1 = EmiAxis(Size-NPtsExX1+1:Size);
    ExFlagX1 = 1;
else
    for j = 1 : NPtsExX1
        ExXDiagSlice1(j) = AbsMData1(Size-NPtsExX1+j,j);
        ExRealXDiagSlice1(j) = RealMData1(Size-NPtsExX1+j,j);
    end
    ExXDiagF1 = EmiAxis(1:NPtsExX1);
    ExFlagX1 = 0;
end

%% Fit slices
% Define starting guess
AEx1 = 400;
dEx = 0.38*meV2Hz;
GEx = 0.15*meV2Hz;
EIS = 0.035*meV2Hz;
EID = 0.020*meV2Hz;
% frac = 1;           % m12/(sqrt(2)*m01) PSF term
param0 = [AEx1 dEx GEx XCor1 YCor1 EIS EID];

% NDipoleMoment = 11;
NParam = 7;         % Include one for residue square
    
param = param0;
dlmwrite(strcat(OutPath,'parameters.dat'),param,'\t');

% Check initial guess
%         PlotSimSpec = 1;
Fit1 = AO_CoCir_Chk_Init_Guess(param,EmiAxis);
Fit1a = Fit1(1:Size);
Fit1 = Fit1(Size+1:end);
Fit1b = Fit1(1:Size);
Fit1 = Fit1(Size+1:end);
Fit1c = Fit1(1:NPtsExD1);
Fit1 = Fit1(NPtsExD1+1:end);
Fit1d = Fit1(1:NPtsExX1);
Fit1 = Fit1(NPtsExX1+1:end);
Fit1e = Fit1(1:Size);
Fit1 = Fit1(Size+1:end);
Fit1f = Fit1(1:Size);
Fit1 = Fit1(Size+1:end);
Fit1g = Fit1;

Fig6 = figure(6);
set(Fig6,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
subplot(2,4,1);
plot(EmiAxis,AbsHSlice1,EmiAxis,Fit1a,'-g');
title('Co-cir Hor');
subplot(2,4,2);
plot(AbsAxis,AbsVSlice1,AbsAxis,Fit1b,'-g');
title('Co-cir Ver');
subplot(2,4,3);
plot(ExDiagF1,ExDiagSlice1,ExDiagF1,Fit1c,'-g');
title('Co-cir Diag');
subplot(2,4,4);
plot(ExXDiagF1,ExXDiagSlice1,ExXDiagF1,Fit1d,'-g');
title('Co-Cir X-diag');
subplot(2,4,5);
plot(EmiAxis,RealHSlice1,EmiAxis,Fit1e,'-g');
title('Co-cir Hor:Real');
subplot(2,4,6);
plot(AbsAxis,RealVSlice1,AbsAxis,Fit1f,'-g');
title('Co-cir Ver:Real');
subplot(2,4,7);
plot(ExXDiagF1,ExRealXDiagSlice1,ExXDiagF1,Fit1g,'-g');
title('Co-Cir X-diag:Real');
saveas(gcf,strcat(OutPath,'Fits0'),'emf');

ContButton = questdlg('Continue?');
if strcmp(ContButton,'Yes')
    Fit = 1;
else
    Fit = 0;
    return;
end

if Fit
    % Fit co-cir Diag(Abs), X-Diag(Abs) and Hor(real) slices
    YData = cat(2,ExDiagSlice1,ExXDiagSlice1,RealHSlice1);
    XData = cat(2,ExDiagF1,ExXDiagF1,EmiAxis);
    ToFit = 1;
    param2 = [param(1) param(2) param(3) param(6) param(7)];
    [param2F,r2,J2] = nlinfit(XData,YData,@AO_MBE_Fit,param2);
    param(1) = param2F(1);
    param(2) = param2F(2);
    param(3) = param2F(3);
    param(6) = param2F(4);
    param(7) = param2F(5);
    dlmwrite(strcat(OutPath,'parameters.dat'),param,'\t');
    ToFit = 0;
    Fit2 = AO_MBE_Fit(param,EmiAxis);
    Fit2a = Fit2(1:NPtsExD1);
    Fit2 = Fit2(NPtsExD1+1:end);
    Fit2b = Fit2(1:NPtsExX1);
    Fit2 = Fit2(NPtsExX1+1:end);
    Fit2c = Fit2;
    Fig7 = figure(7);
    set(Fig7,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    subplot(1,3,1);
    plot(ExDiagF1,ExDiagSlice1,ExDiagF1,Fit2a,'-g');
    title('Co-cir Diag:Abs');
    axis 'square';
    subplot(1,3,2);
    plot(ExXDiagF1,ExXDiagSlice1,ExXDiagF1,Fit2b,'-g');
    axis 'square';
    title('Co-cir X-Diag:Abs');
    subplot(1,3,3);
    plot(EmiAxis,RealHSlice1,EmiAxis,Fit2c,'-g');
    axis 'square';
    title('Co-cir Hor:Real');
    saveas(gcf,strcat(OutPath,'Fits1'),'emf');
    
    % Fit co-cir:all slices starting with best guess from above
    YData = cat(2,AbsHSlice1,AbsVSlice1',ExDiagSlice1,ExXDiagSlice1,...
        RealHSlice1,RealVSlice1',ExRealXDiagSlice1);
    XData = cat(2,EmiAxis,AbsAxis,ExDiagF1,ExXDiagF1,...
        EmiAxis,AbsAxis,ExXDiagF1);
    ToFit = 1;
    param4 = [param(1) param(2) param(3) param(6) param(7)];
    [param4F,r4,J4] = nlinfit(XData,YData,@AO_All_Slices_Fit,param4);
    param(1) = param4F(1);
    param(2) = param4F(2);
    param(3) = param4F(3);
    param(6) = param4F(4);
    param(7) = param4F(5);
    ToFit = 0;
end

%% Final fits
PlotSimSpec = 1;
Fit8 = AO_CoCir_Chk_Init_Guess(param,EmiAxis);
%%
% Normalize Slices
Fit8 = Fit8./VMax1;
AbsHSlice1 = AbsHSlice1./VMax1;
AbsVSlice1 = AbsVSlice1./VMax1;
ExDiagSlice1 = ExDiagSlice1./VMax1;
ExXDiagSlice1 = ExXDiagSlice1./VMax1;
RealHSlice1 = RealHSlice1./VMax1;
RealVSlice1 = RealVSlice1./VMax1;
ExRealXDiagSlice1 = ExRealXDiagSlice1./VMax1;

Fit8a = Fit8(1:Size);
Fit8 = Fit8(Size+1:end);
Fit8b = Fit8(1:Size);
Fit8 = Fit8(Size+1:end);
Fit8c = Fit8(1:NPtsExD1);
Fit8 = Fit8(NPtsExD1+1:end);
Fit8d = Fit8(1:NPtsExX1);
Fit8 = Fit8(NPtsExX1+1:end);
Fit8e = Fit8(1:Size);
Fit8 = Fit8(Size+1:end);
Fit8f = Fit8(1:Size);
Fit8 = Fit8(Size+1:end);
Fit8g = Fit8;

%%
Fig12 = figure(12);
FS = 3.7;       % Fig Size in cm
set(Fig12,'Units', 'Centimeters', 'OuterPosition', [1 1 17.2 15]);
subplot(2,4,1);
plot(ExDiagF1,ExDiagSlice1,'-k',ExDiagF1,Fit8c,'--r','LineWidth',2);
set(gca,'Units','Centimeters','Position',[1 7 FS FS]);
axis([1540 1548 0 1.05]);
set(gca,'YTick',[0 0.5 1],'YTickLabel',[0 0.5 1],...
    'XTick',[1541 1544 1547],'XTickLabel',[1541 1544 1547]);
set(gca,'LineWidth',1,'TickLength',[0.025 0.025]);
title('Co-cir Diag');
subplot(2,4,2);
plot(ExXDiagF1,ExXDiagSlice1,'-k',ExXDiagF1,Fit8d,'--r','LineWidth',2);
set(gca,'Units','Centimeters','Position',[4.8 7 FS FS]);
axis([1540 1548 0 1.05]);
set(gca,'YTick',[0 0.5 1],'YTickLabel',[],...
    'XTick',[1541 1544 1547],'XTickLabel',[1541 1544 1547]);
set(gca,'LineWidth',1,'TickLength',[0.025 0.025]);
title('Co-Cir X-diag');
subplot(2,4,3);
plot(EmiAxis,AbsHSlice1,'-k',EmiAxis,Fit8a,'--r','LineWidth',2);
set(gca,'Units','Centimeters','Position',[8.6 7 FS FS]);
axis([1540 1548 0 1.05]);
set(gca,'YTick',[0 0.5 1],'YTickLabel',[],...
    'XTick',[1541 1544 1547],'XTickLabel',[1541 1544 1547]);
set(gca,'LineWidth',1,'TickLength',[0.025 0.025]);
title('Co-cir Hor');
subplot(2,4,4);
plot(AbsAxis,AbsVSlice1,'-k',AbsAxis,Fit8b,'--r','LineWidth',2);
set(gca,'Units','Centimeters','Position',[12.4 7 FS FS]);
axis([-1548 -1540 0 1.05]);
set(gca,'YTick',[0 0.5 1],'YTickLabel',[],...
    'XTick',[-1547 -1544 -1541],'XTickLabel',[-1547 -1544 -1541]);
set(gca,'LineWidth',1,'TickLength',[0.025 0.025]);
title('Co-cir Ver');
subplot(2,4,6);
plot(ExXDiagF1,ExRealXDiagSlice1,'-k',ExXDiagF1,Fit8g,'--r','LineWidth',2);
set(gca,'Units','Centimeters','Position',[4.8 1 FS FS]);
axis([1540 1548 -0.6 0.8]);
set(gca,'XTick',[1541 1544 1547],'XTickLabel',[1541 1544 1547]);
set(gca,'LineWidth',1,'TickLength',[0.025 0.025]);
title('Co-Cir X-diag:Real');
subplot(2,4,7);
plot(EmiAxis,RealHSlice1,'-k',EmiAxis,Fit8e,'--r','LineWidth',2);
set(gca,'Units','Centimeters','Position',[8.6 1 FS FS]);
axis([1540 1548 -0.6 0.8]);
set(gca,'YTick',[-0.5 0 0.5],'YTickLabel',[],...
    'XTick',[1541 1544 1547],'XTickLabel',[1541 1544 1547]);
set(gca,'LineWidth',1,'TickLength',[0.025 0.025]);
title('Co-cir Hor:Real');
subplot(2,4,8);
plot(AbsAxis,RealVSlice1,'-k',AbsAxis,Fit8f,'--r','LineWidth',2);
set(gca,'Units','Centimeters','Position',[12.4 1 FS FS]);
axis([-1548 -1540 -0.6 0.8]);
set(gca,'YTick',[-0.5 0 0.5],'YTickLabel',[],...
    'XTick',[-1547 -1544 -1541],'XTickLabel',[-1547 -1544 -1541]);
set(gca,'LineWidth',1,'TickLength',[0.025 0.025]);
title('Co-cir Ver:Real');
% saveas(gcf,strcat(OutPath,'FinalFits'),'emf');
saveas(gcf,strcat(OutPath,'FinalFits'),'eps');

%% Convert fit parameter units
param(2) = param(2)./meV2Hz;
param(3) = param(3)./meV2Hz;
param(4) = param(4)./meV2Hz;
param(5) = param(5)./meV2Hz;
param(6) = param(6)./meV2Hz;
param(7) = param(7)./meV2Hz;
dlmwrite(strcat(OutPath,'CoCir_Parameters.dat'),param,...
    'delimiter','\t','precision',6);