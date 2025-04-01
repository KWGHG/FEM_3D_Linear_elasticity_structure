clear all; close all; clc;

figure; hold on; box on; grid on;
view(120,30);

xlim([0,2]);
ylim([0,2]);
zlim([0,2]);

x = [0,1,1,0,0,1,1,0];
y = [0,0,1,1,0,0,1,1];
z = [0,0,0,0,1,1,1,1];
c = [0,0,0,0,1,1,1,1];

patch_queue = [1,2,3,4;5,8,7,6;1,5,6,2;2,6,7,3;3,7,8,4;4,8,5,1];

for i = 1:size(patch_queue,1)
    fill3(x(patch_queue(i,:)),y(patch_queue(i,:)),z(patch_queue(i,:)),...
        c(patch_queue(i,:)));
end

colorbar;
colormap jet;