%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this code aims at providing demonstration of cumulative B_splines for R3 
% SO(3) and SE(3) sets of data 
% author : Mahmoud Z. KHAIRALLAH
% mail: mahmoud.kairallah@univ-evry.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear all; close all;

% read data file
data = readtable('data_gt.csv');
% generate 6-DoF points (R3 and SO(3))
n = 3 ;% (n-1) degree polynomial, namly cubic
I = 100; % number of points
x = table2array(data(1:I,2:4))';
q = table2array(data(1:I,5:8));
q = quatnormalize(q);
p = table2array(data(1:I,2:8));
t = [zeros(1,n-1),linspace(0,1,I-n+2),ones(1,n-1)]; % create uniform knots = n+i
% t = linspace(0,1,size(x,2)+n);
T = linspace(0,1,I);
m1 = 60*I;
m2 = 60;
X = R3_BSplines(n,x,t,m1);


% plot the standard and cumulative spline bases
figure
hold on 
for i = 0:length(t)-n-1
    subplot(2,1,1)
    hold on 
    [YY1,XX1] = spline_basis(i,n,t);
    plot(XX1,YY1)
end
for i = 1:length(t)-2*n+1
    subplot(2,1,2)
    hold on 
    [YY2,XX2] = cumul_spline_basis(i,n,t);
    plot(XX2,YY2)
end



% data using de-boor's formula
[X1,T1] = R3_cumul_BSpline(n,x,t,m1);
[Q1,T2] = SO3_cumul_BSpline(n,q,t,m1);
[P1,Tr,T3] = SE3_cumul_BSpline(n,p,t,m1);

% datausing matrix representation
[X2,T4] = R3_cumul_BSpline_mat(n,x,m2);
[Q2,T5] = SO3_cumul_BSpline_mat(n,q,m2);
[P2,T6] = SE3_cumul_BSpline_mat(n,p,m2);



% showing the data using SE3 splines
tp = theaterPlot('XLimit',[min(x(1,:)) max(x(1,:))],'YLimit',[min(x(2,:)) max(x(2,:))],'ZLimit',[min(x(3,:)) max(x(3,:))]);
op = orientationPlotter(tp,'DisplayName','Splined data',...
    'LocalAxesLength',0.005,'MarkerSize',3,'HistoryDepth',99);
pose_t = quaternion(q);
pose_s = quaternion(P1(:,4:7));
plotOrientation(op,pose_t,x')
pause
for i = 1:length(pose_s)
    plotOrientation(op,pose_s(i),P1(i,1:3))
    drawnow;
end
plotOrientation(op,pose_s(1:100:end),X1(:,1:100:end)')


% plot the 3d space data
figure
plot3(x(1,:),x(2,:),x(3,:),'k.','MarkerSize',10)
hold on 
plot3(X(1,:),X(2,:),X(3,:),'b')
plot3(X1(1,:),X1(2,:),X1(3,:),'r')
plot3(P(:,1),P(:,2),P(:,3),'g')
label = cellstr(num2str([1:I]'));
text(x(1,:),x(2,:),x(3,:),label,'VerticalAlignment','bottom','HorizontalAlignment','right')
grid on


% plot each component of the quaternions alone
figure
subplot(4,1,1)
hold on 
plot(T,q(:,1))
plot(T2,Q(:,1))
plot(T3,P(:,4))
subplot(4,1,2)
hold on 
plot(T,q(:,2))
plot(T2,Q(:,2))
plot(T3,P(:,5))
subplot(4,1,3)
hold on 
plot(T,q(:,3))
plot(T2,Q(:,3))
plot(T3,P(:,6))
subplot(4,1,4)
hold on 
plot(T,q(:,4))
plot(T2,Q(:,4))
plot(T3,P(:,7))


