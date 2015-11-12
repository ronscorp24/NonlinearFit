function F = AO_CoCir_Chk_Init_Guess(param,~)

global AbsFMesh EmiFMesh
global tMesh tauMesh
global ExPkX1 ExPkY1
global NPtsExD1 NPtsExX1 ExFlagD1 ExFlagX1
global Size
global PlotSimSpec OutPath
% global meV2Hz

AEx1 = param(1);
dEx = param(2);
GEx = param(3);
XCor1 = param(4);
YCor1 = param(5);
EIS = param(6);
EID = param(7);
frac = 1;

% Co-circular spectrum
Sig_a = exp(-(GEx.*(tMesh+tauMesh))).*...
    exp(-(1/2)*dEx^2.*(tMesh - tauMesh).^2);
% Norm_a = sqrt(sum(sum(abs(Sig_a).^2)));
Sig_b = exp(1i*EIS.*tMesh).*...
    exp(-(GEx.*tauMesh + (GEx+EID).*tMesh)).*...
    exp(-(1/2)*dEx^2.*(tMesh - tauMesh).^2);
% Norm_b = sqrt(sum(sum(abs(Sig_b).^2)));
S1SigTime1 = 2*exp(1i*XCor1.*tMesh).*exp(1i*YCor1.*tauMesh).*AEx1.*...
    (Sig_a - frac.^2.*Sig_b);
S1SigOmega1 = fftshift(fft(fftshift(fft(S1SigTime1,[],2),2),[],1),1);
S1AbsSpec1 = abs(S1SigOmega1);
S1RealSpec1 = real(S1SigOmega1);

if PlotSimSpec
    Max1 = max(max(S1AbsSpec1));
    Fig5 = figure(5);
    set(Fig5,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    subplot(121);
    contourf(EmiFMesh, AbsFMesh, S1AbsSpec1, linspace(0,Max1,50),...
        'LineStyle', 'none');
    line([EmiFMesh(1) EmiFMesh(end)], [AbsFMesh(end) AbsFMesh(1)],...
       'LineStyle', '--', 'Color', [1 1 1],'LineWidth',1.5);
    axis square;
    subplot(122)
    contourf(EmiFMesh, AbsFMesh, S1RealSpec1, linspace(-Max1,Max1,50),...
        'LineStyle', 'none');
    line([EmiFMesh(1) EmiFMesh(end)], [AbsFMesh(end) AbsFMesh(1)],...
       'LineStyle', '--', 'Color', [1 1 1],'LineWidth',1.5);
    axis square;
    saveas(gcf,strcat(OutPath,'FitSpectra'),'emf');
    
    % Save fit 2D spectra data
    dlmwrite(strcat(OutPath,'CoCirAbsSpecFit.dat'),S1AbsSpec1,'Delimiter','\t','Precision',6);
    dlmwrite(strcat(OutPath,'CoCirRealSpecFit.dat'),S1RealSpec1,'Delimiter','\t','Precision',6);
end

% Horizontal Slices
% Co-circular
AbsHSliceFit1 = S1AbsSpec1(ExPkY1,:);
RealHSliceFit1 = S1RealSpec1(ExPkY1,:);

% Vertical slices
% Co-circular
AbsVSliceFit1 = S1AbsSpec1(:,ExPkX1);
RealVSliceFit1 = S1RealSpec1(:,ExPkX1);

% Diagonal slices
% Co-circular
AbsDSliceFit1 = zeros(1,NPtsExD1);
if ExFlagD1
    for j = 1:NPtsExD1
        AbsDSliceFit1(j) = S1AbsSpec1(NPtsExD1+1-j,j);
    end
else
    for j = 1:NPtsExD1
        AbsDSliceFit1(j) = S1AbsSpec1(Size+1-j,Size-NPtsExD1+j);
    end
end

% X-Diagonal slices
% Co-circular
AbsXDSliceFit1 = zeros(1,NPtsExX1);
RealXDSliceFit1 = zeros(1,NPtsExX1);
if ExFlagX1
    for j = 1 : NPtsExX1
        AbsXDSliceFit1(j) = S1AbsSpec1(j,Size-NPtsExX1+j);
        RealXDSliceFit1(j) = S1RealSpec1(j,Size-NPtsExX1+j);
    end
else
    for j = 1 : NPtsExX1
        AbsXDSliceFit1(j) = S1AbsSpec1(Size-NPtsExX1+j,j);
        RealXDSliceFit1(j) = S1RealSpec1(Size-NPtsExX1+j,j);
    end
end

F = cat(2,AbsHSliceFit1,AbsVSliceFit1',AbsDSliceFit1,AbsXDSliceFit1,...
    RealHSliceFit1,RealVSliceFit1',RealXDSliceFit1);