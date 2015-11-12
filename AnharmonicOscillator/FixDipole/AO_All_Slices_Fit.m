function F = AO_All_Slices_Fit(p0,~)

global OutPath
global tMesh tauMesh
global ExPkX1 ExPkY1
global NPtsExD1 NPtsExX1 ExFlagD1 ExFlagX1
global Size
global ToFit
param = dlmread(strcat(OutPath,'parameters.dat'),'\t');

AEx1 = param(1);
dEx = param(2);
GEx = param(3);
XCor1 = param(4);
YCor1 = param(5);
EIS = param(6);
EID = param(7);
frac = 1;

if ToFit
    AEx1 = p0(1);
    dEx = p0(2);
    GEx = p0(3);
    EIS = p0(4);
    EID = p0(5);
%     frac = p0(6);
end

% Co-circular spectrum
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