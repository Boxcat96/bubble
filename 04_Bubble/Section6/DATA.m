clc;
close all;

if pick_reload == 1
    data = readmatrix('..\rawdata\HP.csv');
    data = data(2:end, 2:end)';

    zz = data(1,:);
    zzcredit = data(2,:);
    zzstock = data(3,:);

    save DATA
end

load DATA



