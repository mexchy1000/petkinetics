function DVRmap=logan_map(imgdata,ref,midtime,startframe,endframe)

%imgdata: 4D-matrix (x,y,z,t)
% ref: time x counts values of reference ROI (frames x 1 vector)
% midtime: Time frame value.
% startframe, endframe: parameters for which frames involved
% midtime : min

zsize=size(imgdata,3);
DVRmap=ones(size(imgdata,1),size(imgdata,2),size(imgdata,3));
Thrs=0.5; % Threshold for calcuation. ratio of voxel counts to sum of reference counts 

for zval=1:zsize
    disp('slice no :');
    disp(zval);
    for i=1:size(imgdata,1)
        for j=1:size(imgdata,2)
            data=imgdata(i,j,zval,:);
            data=squeeze(data);
            if sum(data)< 0.5*sum(ref)
                DVR(i,j,zval)=1;
            else
                DVR(i,j,zval)=logan_single(data,ref,midtime,startframe,endframe);
            end
        end
    end
end

DVRmap=DVR;

            
            
    