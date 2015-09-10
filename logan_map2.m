function DVRmap=logan_map2(imgdata,ref,midtime,startframe,endframe)

%imgdata: 4D-matrix (x,y,z,t)
% ref: time x counts values of reference ROI (frames x 1 vector)
% midtime: Time frame value.
% startframe, endframe: parameters for which frames involved
% midtime : min
% 20150909 by Choi H.

zsize=size(imgdata,3);
DVRmap=ones(size(imgdata,1),size(imgdata,2),size(imgdata,3));
Thrs=0.5; % Threshold for calcuation. ratio of voxel counts to sum of reference counts 

for zval=1:zsize
    disp('slice no :');
    disp(zval); % slice no check
    
    data=imgdata(:,:,zval,:);
    data=reshape(data,[],size(imgdata,4)); % matrix form : voxel counts x frames
    data=data';
    
    X=zeros(endframe,size(data,2));
    Y=zeros(endframe,size(data,2));
    DVRvec=ones(size(data,2),1);
    
    refmat=repmat(ref,[1 size(data,2)]);
    
    X(1,:)=0.5*midtime(1).*refmat(1,:)./data(1,:);
    Y(1,:)=0.5*midtime(1).*data(1,:)./data(1,:);
    
    for i=2:endframe
        X(i,:)=trapz(midtime(1:i),refmat(1:i,:))./data(i,:);
        Y(i,:)=trapz(midtime(1:i),data(1:i,:))./data(i,:);
    end
    
    X=X([startframe:endframe],:);
    Y=Y([startframe:endframe],:);
    
    %estimation of coefficient
    for i=1:size(data,2)
        if sum(data(:,i))>0.5*sum(refmat(:,i))
            Xvec=X(:,i);
            Yvec=Y(:,i);
            coeff=[ones(length(Xvec),1) Xvec]\Yvec;
            DVRvec(i)=coeff(2);
        end
    end
    
    DVR=reshape(DVRvec,size(imgdata,1),size(imgdata,2));
    DVRmap(:,:,zval)=DVR;
end

    

            
            
    