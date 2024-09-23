load("ELE532_Lab1_Data.mat")
y_audio = x_audio;
rows = size(y_audio,1);
cols = size(y_audio,2);
threshold = 0.01;
num_of_zeros = 0;
for i = 1:rows
    for j = 1:cols
        if(abs(y_audio(i,j))<threshold)
            y_audio(i,j) = 0;
            num_of_zeros = num_of_zeros + 1;
        end
    end
end
fprintf("\nNumber of zeros: " + num_of_zeros +"\n");
sound(y_audio,8000);