%% cascated image
function [O,T_Average,NEW] = main_V3_function_patch(path,Pic_num,ori)

[row, col,channel] = size(ori);
O = zeros(row,col,channel);
T_Average = zeros(row,col,channel);
NEW = zeros(row,col,channel);

for k = 1:channel

    sigma = 7;
    len = 21;  %odd number;
    omega = 2*pi/(len+1);
    modify_cos_line = cos(omega*(0:len-1)).*exp(-((0:len-1)-len/2).^2/sigma.^2);
    modify_sin_line = sin(omega*(0:len-1)).*exp(-((0:len-1)-len/2).^2/sigma.^2);

    modify_cos_patch = zeros(len,len);
    modify_sin_patch =  zeros(len,len);

    for i = 1:len
        modify_cos_patch(i,:) = modify_cos_line;
        modify_sin_patch(i,:) = modify_sin_line;
    end

    NV = zeros(row,col,Pic_num-1);
    Diff = zeros(row,col,Pic_num-1);

    for z = 1:Pic_num-1 %1:Pic_num-1

        z

        name = [int2str(z),'.jpg'];%.jpg png
        L = double(imread([path,name]));
        L = L(:,:,k);

        name = [int2str(z+1),'.jpg'];%.jpg
        R = double(imread([path,name]));
        R = R(:,:,k);

        clear name;

        for i = 1:row
            for j  = 1:col

                if(i<floor(len/2)+1)

                    if(j<floor(len/2)+1)   %max_j = 25
                        I1_sin = L(1:i+floor(len/2),1:j+floor(len/2)).*modify_sin_patch(1:floor(len/2)+i,floor(len/2)+2-j:2*floor(len/2)+1);     %floor(len/2)+j
                        I1_cos = L(1:i+floor(len/2),1:j+floor(len/2)).*modify_cos_patch(1:floor(len/2)+i,floor(len/2)+2-j:2*floor(len/2)+1);   %floor(len/2)+j
                        I2_sin = R(1:i+floor(len/2),1:j+floor(len/2)).*modify_sin_patch(1:floor(len/2)+i,floor(len/2)+2-j:2*floor(len/2)+1);
                        I2_cos = R(1:i+floor(len/2),1:j+floor(len/2)).*modify_cos_patch(1:floor(len/2)+i,floor(len/2)+2-j:2*floor(len/2)+1);
                    elseif (j>col-floor(len/2))
                        I1_sin = L(1:i+floor(len/2), j-floor(len/2):col).*modify_sin_patch(1:i+floor(len/2),1:floor(len/2)+col-j+1 );
                        I1_cos = L(1:i+floor(len/2),j-floor(len/2):col).*modify_cos_patch(1:i+floor(len/2),1:floor(len/2)+col-j+1 );
                        I2_sin = R(1:i+floor(len/2), j-floor(len/2):col).*modify_sin_patch(1:i+floor(len/2),1:floor(len/2)+col-j+1 );
                        I2_cos = R(1:i+floor(len/2),j-floor(len/2):col).*modify_cos_patch(1:i+floor(len/2),1:floor(len/2)+col-j+1 );
                    else
                        I1_sin = L(1:i+floor(len/2),j-floor(len/2):j+floor(len/2)).*modify_sin_patch(1:i+floor(len/2),:);
                        I1_cos = L(1:i+floor(len/2),j-floor(len/2):j+floor(len/2)).*modify_cos_patch(1:i+floor(len/2),:);
                        I2_sin = R(1:i+floor(len/2),j-floor(len/2):j+floor(len/2)).*modify_sin_patch(1:i+floor(len/2),:);
                        I2_cos = R(1:i+floor(len/2),j-floor(len/2):j+floor(len/2)).*modify_cos_patch(1:i+floor(len/2),:);
                    end

                elseif(i>row-floor(len/2))

                    if(j<floor(len/2)+1)   %max_j = 25
                        I1_sin = L(1:floor(len/2)+row-i+1,1:j+floor(len/2)).*modify_sin_patch(1:floor(len/2)+row-i+1,floor(len/2)+2-j:2*floor(len/2)+1);     %floor(len/2)+j
                        I1_cos = L(1:floor(len/2)+row-i+1,1:j+floor(len/2)).*modify_cos_patch(1:floor(len/2)+row-i+1,floor(len/2)+2-j:2*floor(len/2)+1);   %floor(len/2)+j
                        I2_sin = R(1:floor(len/2)+row-i+1,1:j+floor(len/2)).*modify_sin_patch(1:floor(len/2)+row-i+1,floor(len/2)+2-j:2*floor(len/2)+1);
                        I2_cos = R(1:floor(len/2)+row-i+1,1:j+floor(len/2)).*modify_cos_patch(1:floor(len/2)+row-i+1,floor(len/2)+2-j:2*floor(len/2)+1);
                    elseif (j>col-floor(len/2))
                        I1_sin = L(1:floor(len/2)+row-i+1, j-floor(len/2):col).*modify_sin_patch(1:floor(len/2)+row-i+1,1:floor(len/2)+col-j+1 );
                        I1_cos = L(1:floor(len/2)+row-i+1,j-floor(len/2):col).*modify_cos_patch(1:floor(len/2)+row-i+1,1:floor(len/2)+col-j+1 );
                        I2_sin = R(1:floor(len/2)+row-i+1, j-floor(len/2):col).*modify_sin_patch(1:floor(len/2)+row-i+1,1:floor(len/2)+col-j+1 );
                        I2_cos = R(1:floor(len/2)+row-i+1,j-floor(len/2):col).*modify_cos_patch(1:floor(len/2)+row-i+1,1:floor(len/2)+col-j+1 );
                    else
                        I1_sin = L(1:floor(len/2)+row-i+1,j-floor(len/2):j+floor(len/2)).*modify_sin_patch(1:floor(len/2)+row-i+1,:);
                        I1_cos = L(1:floor(len/2)+row-i+1,j-floor(len/2):j+floor(len/2)).*modify_cos_patch(1:floor(len/2)+row-i+1,:);
                        I2_sin = R(1:floor(len/2)+row-i+1,j-floor(len/2):j+floor(len/2)).*modify_sin_patch(1:floor(len/2)+row-i+1,:);
                        I2_cos = R(1:floor(len/2)+row-i+1,j-floor(len/2):j+floor(len/2)).*modify_cos_patch(1:floor(len/2)+row-i+1,:);
                    end

                else

                    if(j<floor(len/2)+1)   %max_j = 25
                        I1_sin = L(i-floor(len/2):i+floor(len/2),1:j+floor(len/2)).*modify_sin_patch(:,floor(len/2)+2-j:2*floor(len/2)+1);     %floor(len/2)+j
                        I1_cos = L(i-floor(len/2):i+floor(len/2),1:j+floor(len/2)).*modify_cos_patch(:,floor(len/2)+2-j:2*floor(len/2)+1);   %floor(len/2)+j
                        I2_sin = R(i-floor(len/2):i+floor(len/2),1:j+floor(len/2)).*modify_sin_patch(:,floor(len/2)+2-j:2*floor(len/2)+1);
                        I2_cos = R(i-floor(len/2):i+floor(len/2),1:j+floor(len/2)).*modify_cos_patch(:,floor(len/2)+2-j:2*floor(len/2)+1);
                    elseif (j>col-floor(len/2))
                        I1_sin = L(i-floor(len/2):i+floor(len/2), j-floor(len/2):col).*modify_sin_patch(:,1:floor(len/2)+col-j+1 );
                        I1_cos = L(i-floor(len/2):i+floor(len/2),j-floor(len/2):col).*modify_cos_patch(:,1:floor(len/2)+col-j+1 );
                        I2_sin = R(i-floor(len/2):i+floor(len/2), j-floor(len/2):col).*modify_sin_patch(:,1:floor(len/2)+col-j+1 );
                        I2_cos = R(i-floor(len/2):i+floor(len/2),j-floor(len/2):col).*modify_cos_patch(:,1:floor(len/2)+col-j+1 );
                    else
                        I1_sin = L(i-floor(len/2):i+floor(len/2),j-floor(len/2):j+floor(len/2)).*modify_sin_patch(:,:);
                        I1_cos = L(i-floor(len/2):i+floor(len/2),j-floor(len/2):j+floor(len/2)).*modify_cos_patch(:,:);
                        I2_sin = R(i-floor(len/2):i+floor(len/2),j-floor(len/2):j+floor(len/2)).*modify_sin_patch(:,:);
                        I2_cos = R(i-floor(len/2):i+floor(len/2),j-floor(len/2):j+floor(len/2)).*modify_cos_patch(:,:);
                    end

                end

                I1_sinfft = fft2(I1_sin);
                I2_sinfft = fft2(I2_sin);
                I1_cosfft = fft2(I1_cos);
                I2_cosfft = fft2(I2_cos);

                atanI_Sin = atan2(imag(I1_sinfft-I2_sinfft),real(I1_sinfft-I2_sinfft));
                atanI_Cos = atan2(imag(I1_cosfft-I2_cosfft),real(I1_cos-I2_cosfft));

                theta = mean(mean(atanI_Cos))/2+mean(mean(atanI_Sin))/2;

                %NewH(i,j) = sqrt( 1*L(i,j).^2 + 1*R(i,j).^2+2*cos(mean(theta)).*(L(i,j)*R(i,j)) )/sqrt(4);
                NewH_cos(i,j) = sqrt(abs(cos(mean(theta))).*(L(i,j)*R(i,j)));


            end
        end



        %aug0 = find( abs(NewV_cos/2+NewH_cos/2-(L/2+R/2))>1 );
        %NewV_cos(aug0) = 0;
        NV(:,:,z) = NewH_cos;
        Diff(:,:,z)  = NewH_cos-L/2-R/2;
    end


    T_average = zeros(row,col);
    for i = 1:Pic_num %1:Pic_num
        name = [int2str(i),'.jpg'];%.jpg
        I = double(imread([path,name]));
        I = I(:,:,k);  
        T_average = T_average + I;
    end
    T_average  = T_average/Pic_num;

    Diff1 = Diff;
    NV1 = NV;
    for i = 1:Pic_num-1 %1:Pic_num-1
        aug0 = Diff1(:,:,i);
        aug1 = find(abs(aug0)>1);
        aug2 = NV1(:,:,i);
        aug2(aug1) = 0;
        NV1(:,:,i) = aug2;
    end

    New = zeros(row,col);
    for i = 1:row
        for j = 1:col
            aug0 = find(NV1(i,j,:)==0);
            if length(aug0) == Pic_num-1
                continue;
            end
            coff = Pic_num-1-length(aug0);
            New(i,j) = sum(NV1(i,j,:))/coff;
        end
    end
    aug0 = find(New ==0 );
    New(aug0)  = T_average(aug0);

    o = New-T_average;
    
    O(:,:,k) = o;
    NEW(:,:,k) = New;
    T_Average(:,:,k) = T_average;
    
end
