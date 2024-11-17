clc;clear
for ii = 0 : 8

path = ['..\',num2str(ii),'\'];
ori = imread([path,'OriginalPic.png']) ;
% ori = imread([path,'1.jpg']) ;

%% GI is a index to measure the complexity of the scene. 
%% Generally speaking, a good improvement would be obtained if GI>0.45
[GI,~,~] = HoGVariety(ori,18);  
disp(GI);
Pic_num = 40;
[feature,average,New] = main_V3_function_patch(path,Pic_num,ori(:,:,:));

ori = ori(:,:,:);

%figure;imshow(uint8(ori),'Border','Tight')
%mesh((feature(:,:,1)));axis off;view([0 270]);colormap(jet);colorbar
%figure;imshow(uint8(T_average),'Border','Tight')
%figure;imshow((uint8(New+1*feature)),'Border','Tight');

imwrite(ori,[path,'OriginalPic.png']);
imwrite(uint8(average),[path,'average.png']);
imwrite((uint8(New+1*feature)),[path,'BIPA.png']);


fileID = fopen([path,'data.txt'], 'w');

a = 1;
Example = imread([path,'1.jpg']);
NF = New + a * feature;
fprintf('\n')

for i = 1:3
    O = ori(:,:,i);
    nf = NF(:,:,i);
    TA = average(:,:,i);
    e = Example(:,:,i);

    % 写入 SSIM 结果到文件
    fprintf(fileID, 'SSIM between original and NF for channel %d: %0.4f\n', i, ssim(double(O), double(nf)));
    fprintf(fileID, 'SSIM between original and T_average for channel %d: %0.4f\n', i, ssim(double(O), double(TA)));
    fprintf(fileID, 'SSIM between original and Example for channel %d: %0.4f\n', i, ssim(double(O), double(e)));
end

fclose(fileID);

end