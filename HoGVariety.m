function [GI,Angle,Amplititude] = HoGVariety(Img,bin)

[row,col,channel] = size(Img);
vector = zeros(channel,bin);

for k = 1:channel

    img = double(Img(:,:,k))/255.0;

    Amplititude = zeros(row-1,col-1);
    Angle = Amplititude;


    for i = 2:row - 1   %sobel边缘检测
        for j = 2:col - 1

            Gx = (img(i+1,j-1) + 2*img(i+1,j) + img(i+1,j+1)) - (img(i-1,j-1) + 2*img(i-1,j) + img(i-1,j+1)); %Horizontal Direction
            Gy = (img(i-1,j+1) + 2*img(i,j+1) + img(i+1,j+1)) - (img(i-1,j-1) + 2*img(i,j-1) + img(i+1,j-1)); %Vertical Direction
            Amplititude(i,j) = sqrt(Gx^2 + Gy^2);
            Angle(i,j) = abs(atan2(Gy,Gx));

            if Angle(i,j)<0
                Angle(i,j) =  Angle(i,j)+pi;
            end

            bin_num = floor(Angle(i,j)/pi*bin);

            if  bin_num == bin
                vector(k,1) = vector(k,1)+Amplititude(i,j);
            else
                vector(k,bin_num+1) = vector(k,bin_num+1)+Amplititude(i,j);
            end

        end
    end

end
vector = mean(vector,1);
GI = sum(vector/row/col);


