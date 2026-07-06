%compute longitudinal equilibrium points

%% choose drone to trim
drone = generateDarko();

%% definition of pitch angles to test
thetas = (2:2:110)*pi/180;
N = length(thetas);

%% init of polars
vinfs = zeros(1,N);
omegs = zeros(1,N);
dels  = zeros(1,N);
Ns    = zeros(1,N);

%% computation of polars by iteration on alphas
for i=1:N
    [xts, uts, No] = longTrim( thetas(i), drone );
    if No > 0
        vinfs(i) = xts(1);
        omegs(i) = uts(2);
        dels(i)  = uts(3);
        Ns(i)    = No;
    else
        vinfs(i) = 0;
        omegs(i) = 0;
        dels(i)  = 0;
        Ns(i)    = No;
    end
end

%% windtunnel longitudinal polars
% THIS DATA IS NOT COMPLETE? ONLY A SUBSET!!
pol_alpha = fliplr([ 6 7 10 15 20 30 40 50 60 70 80 ]);
pol_vinf = fliplr([ 19.2 15.7 13.5 11.1 10.3 8.2 7.0 4.5 2.8 1.7 0.0 ]);
pol_pe_delta1 = fliplr([ 0.9775  1.9550   2.4438  4.3988  8.7977  18.8172  24.6823  20.0391  12.4633  4.1544  2.1994  ]);
pol_pe_delta2 = fliplr([ 3.1769  1.7107   3.4213  4.6432  12.4633 22.4829  28.8368  23.9492  14.4184  4.3988  2.9326  ]);
pol_delta = -(pol_pe_delta1 + pol_pe_delta2)/2;
pol_omega = fliplr([ 915  755 674 630 670 684 700 732 764 785 820  ]);

%% plot results
% figure;
% subplot(2,2,1);
% plot(thetas*180/pi,vinfs); hold;
% plot(pol_alpha,pol_vinf,'-x');
% subplot(2,2,2);
% plot(thetas*180/pi,omegs); hold;
% plot(pol_alpha,pol_omega,'-x');
% subplot(2,2,3);
% plot(thetas*180/pi,dels*180/pi); hold;
% plot(pol_alpha,pol_delta,'-x');
% subplot(2,2,4);
% plot(thetas*180/pi,Ns);

%% export results to latex
M1 = [180/pi*thetas' vinfs' omegs' 180/pi*dels'];
%save('..\ieeeconf\figs\polaires\polTheory.dat','M1','-ascii');
M2 = [pol_alpha' pol_vinf' pol_omega' pol_delta'];
%save('..\ieeeconf\figs\polaires\polTunnel.dat','M2','-ascii');

%% solve for lateral polars as well!
% list of pilot commands to go through
v0s = 0:1:20; V = length(v0s); % 2 - 0:1:20
w0s = 0:5*pi/180:50*pi/180; W = length(w0s); % & - 0:5*pi/180:50*pi/180
% init data structures for energy-maneuverability diagrams
thets = zeros(V,W);
phis   = zeros(V,W);
w1s    = zeros(V,W);
w2s    = zeros(V,W);
d1s    = zeros(V,W);
d2s    = zeros(V,W);
fvals  = zeros(V,W);
% compute energy maneuverability diagram!
for j = 1:W
    for i = 1:V
        % in case we are in longitudinal mode
        if j == 1
            % we use this value for longitudinal search
            x0 = [45*pi/180; 0; 500; -500; 0; 0];
        else
            % if not we use the one computed before for a slightly
            % less angular velocity
            x0 = [ thets(i,j-1); phis(i,j-1); w1s(i,j-1); w2s(i,j-1); d1s(i,j-1); d2s(i,j-1) ];
        end
        % trim computation
        [x,fval] = fsolve(@(x)latCost(x,v0s(i),w0s(j), drone), x0);
        thets(i,j) = x(1);
        phis(i,j)   = x(2);
        w1s(i,j)    = x(3);
        w2s(i,j)    = x(4);
        d1s(i,j)    = x(5);
        d2s(i,j)    = x(6);
        fvals(i,j)  = min(norm(fval),1);
    end
end

%% filter out divergent trim points
thets_f  = zeros(V,W);
phis_f   = zeros(V,W);
w1s_f    = zeros(V,W);
w2s_f    = zeros(V,W);
d1s_f    = zeros(V,W);
d2s_f    = zeros(V,W);
for j = 1:W
    for i = 1:V
        % if evaluation of function is too big..
        if fvals(i,j) > 0.00001
            thets_f(i,j)  = 10;
            phis_f(i,j)   = 10;
            w1s_f(i,j)    = 1000;
            w2s_f(i,j)    = -1000;
            d1s_f(i,j)    = 10;
            d2s_f(i,j)    = 10;
        else
            thets_f(i,j)  = thets(i,j);
            phis_f(i,j)   = phis(i,j);
            w1s_f(i,j)    = w1s(i,j);
            w2s_f(i,j)    = w2s(i,j);
            d1s_f(i,j)    = d1s(i,j);
            d2s_f(i,j)    = d2s(i,j);
        end
    end
end

%% clean up trim points
for i=1:V
    for j=2:W
        if fvals(i,j) > 0.00001
            thets_f(i,j)  = thets(i,j-1);
            phis_f(i,j)   = phis(i,j-1);
            w1s_f(i,j)    = w1s(i,j-1);
            w2s_f(i,j)    = w2s(i,j-1);
            d1s_f(i,j)    = d1s(i,j-1);
            d2s_f(i,j)    = d2s(i,j-1);
        end
    end
end