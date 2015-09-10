function []=logan_parametricimage(imfile,ref,midtime,startframe,endframe)

% Make parametric image using logan plot
% Imput : imfile - hdr/img format (analyze format, 4D file)
% ref: reference ROI counts or plasma counts
%midtime: minutes
% start/endframe: parameters for which frames involved
% 20150909by Choi H.

fnm=strcat('DVR_',imfile);

imgdata=readanalyze2(imfile);
fhead=analyze75info(imfile);
disp('image files are loaded.');
disp('calculating logan plot...');
DVRmap=logan_map2(imgdata,ref,midtime,startframe,endframe);

writeanalyze2(DVRmap,[size(imgdata,1) size(imgdata,2) size(imgdata,3)],fnm,fhead.PixelDimensions);
disp('parametric images are saved');