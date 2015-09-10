function [DVR,X,Y]=logan_single(data,ref,midtime,startframe,endframe)

% data:  counts values
% ref: time x counts values of reference ROI (frames x 1 vector)
% midtime: Time frame value.
% startframe, endframe: parameters for which frames involved
% dT : min
% Copyright Choi H 09092015


X=zeros(endframe,1);
Y=zeros(endframe,1);

X(1)=0.5*midtime(1)*ref(1)./data(1);
Y(1)=0.5*midtime(1)*data(1)./data(1);

for i=2:endframe
    X(i)=trapz(midtime(1:i),ref(1:i))./data(i);
    Y(i)=trapz(midtime(1:i),data(1:i))./data(i);
end

X2=X(startframe:endframe);
Y2=Y(startframe:endframe);

coeff=[ones(length(X2),1) X2]\Y2;
DVR=coeff(2);
b=coeff(1);

figure;
scatter(X,Y,'.','k');
hold on
plot(X,X*DVR+b,'--');
grid on
xlabel('Int(Cref)/Ct');
ylabel('Int(Ct)/Ct');
title('Logan Plot');

