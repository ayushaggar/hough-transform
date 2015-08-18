% Ayush Aggarwal
% 11177
% Predicting Planes using hough transform
clc;
clear all;
D = textread('CE676_HA7.txt','','delimiter',',' ); 
X=D(:,1);
Y=D(:,2);
Z=D(:,3);
leng=length(D(:,1));

%%
% variables controlling plane search
dist_incr = 0.03;
theta_incr = 2;
phi_incr=2;
maxx=max(D(:, 1));
minx=min(D(:, 1));
maxy=max(D(:, 2));
miny=min(D(:, 2));
maxz=max(D(:, 3));
minz=min(D(:, 3));
Dis_min= sqrt(minx^2+miny^2+minz^2);
Dis_max= sqrt(maxx^2+maxy^2+maxz^2);


%Setting up Houghroom matrix based on size of the D set
% we are taking 3d variable as distance from origin theta and phi
distance = Dis_min: dist_incr : Dis_max;
theta = 0 : theta_incr : 90-theta_incr;
phi= 0:phi_incr:90-phi_incr;   % as x y z are positive
theta_rad = theta*pi/180;
phi_rad=phi*pi/180;
c_phi=cos(phi_rad);
t_theta= length(theta);
t_phi=length(phi);
t_distance=length(distance);
ratio= (t_distance-1)/(distance(t_distance)-distance(1));
Hough_room = zeros(length(theta), length(phi),length(distance)); % making matrix for putting hough values

%%
% if parameters lie with in various parameter then increse the hough room count
for l=1:1:t_distance    
   for k=1:1:t_phi
     for j=1:1:t_theta       
        for i = 1:leng
            dist = X(i)*(cos(phi_rad(k))*cos(theta_rad(j))) + Y(i)*(cos(phi_rad(k)).*sin(theta_rad(j)))+ Z(i)*(sin(phi_rad(k))) ;
            % taking error in distamce = error in point given in question for less
            % computation so hough interval decided by it
            if  distance(l)<(dist+0.015) & distance(l)> (dist-0.015)
                 Hough_room(k, j, l) = Hough_room(k, j, l) +1;
            end
        end
     end
   end
end

%% Plot hough room
plot3(Hough_room(:,1), Hough_room(:,2), Hough_room(:,3));

%%
%For hough room thresolding
thresh = 0.5 *(max(max(Hough_room(:))));

% getting points of result by hough room thershold
[r_theta r_phi r_distance] = find(Hough_room > thresh);

r_Houghroom = Hough_room - thresh;
r_hough_dist = [];
r_hough_theta = [];
r_hough_phi=[];

%Finding x, y, z of the cloud in result Hough_room after thershold.
count=0;
for l = 1:t_distance
    for j=1:t_theta
        for k=1:t_phi
        if r_Houghroom(j, k, l) >= 0
            count=count+1;
            r_hough_dist(count,1) = distance(l);
            r_hough_theta(count,1) = theta(j);
            r_hough_phi(count,1) = phi(k);
        end
        end
    end
end

% Calculation of lines.
r_hough_dist = r_hough_dist * dist_incr;
r_hough_theta = (r_hough_theta * theta_incr) - theta_incr;
r_hough_phi = (r_hough_phi*phi_incr)- phi_incr;
line=[r_hough_dist r_hough_theta r_hough_phi];
dlmwrite('result_line.txt',line,'delimiter',' ','precision','%.6f');
