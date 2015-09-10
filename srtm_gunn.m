function [BP,RI,k2,modfit]=srtm_gunn(scantime,Cexp,Cr,nIter)

%%srtm_gunn by Choi H (2013.8.14)
%%scantime : 2xm matrix, 1st column for frame start (min) and 2nd column
%%for frame end (min)
%%Cexp : Experimental data, receptor rich region
%% Cr : reference data
%% nIter: Iteration numbers. recommend for >1000
%% 

nFrames = size(Cexp, 1);
nTACS = size(Cexp, 2);
beginF = scantime(:,1);
endF = scantime(:,2);
midFrame = (beginF  + endF ) / 2;


frameDur = endF - beginF;


frameCounts = ones([nFrames 1]);
weights = frameDur .^ 2 ./ frameCounts;
W = diag(sqrt(weights));

theta3 = logspace(-4, -1, nIter) * 60; 
B = zeros([nFrames nIter]);
% integral_0^t (Cr(t-tau) * exp(-theta3*tau)) dtau
for i = 1:nIter,
    B(1,i)=0;
    for j = 2:nFrames,
        B(j,i)=trapz(midFrame(1:j),Cr(1:j).*exp(-theta3(i).*(midFrame(j)-midFrame(1:j))));
    end
end

M = zeros([2*nIter nFrames]);
for i = 1:nIter,
    A = [Cr B(:,i)];
    [Q, R] = qr(W * A);
    M(2*(i-1)+[1 2], :) = R \ (Q.');
end

RI = zeros([1 nTACS]);
k2 = zeros([1 nTACS]);
BP = zeros([1 nTACS]);
Ct = zeros([nFrames nIter]);
theta = zeros([2 nIter]);
 modfit.ModelY = zeros([nFrames nTACS]);
 for j = 1:nTACS,
    CPET = Cexp(:, j);
for i = 1:nIter,
    theta(:,i) = M(2*(i-1)+[1 2], :) * W * CPET;
    Ct(:,i) = theta(1,i) * Cr + theta(2,i) * B(:,i);
end
RSS = sum(repmat(weights, [1 nIter]) .* ((repmat(CPET, [1 nIter]) - Ct) .^ 2), 1);
% semilogx(theta3,RSS); xlabel('Theta 3'); ylabel('RSS');
[RSSmin, imin] = min(RSS);
RI(j) = theta(1, imin);
k2(j) = theta(2, imin) + RI(j)*(theta3(imin));
BP(j) = k2(j) / (theta3(imin)) - 1;
modfit.ModelY(:, j) = Ct(:,imin)
 end

end
