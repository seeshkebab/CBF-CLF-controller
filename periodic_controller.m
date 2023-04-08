close all;
clear all;
%% parameters

dt = 0.02; %time step
sim_t = 30; %simulation time period

cbf_gamma0 = 1;
L = 1;

% Initial state.
x0 = [6; 5];

all_roots_vector_secant = [];
all_roots_vector_Bisection = [];

x1d = -7;
x2d = 0;


%% define system

syms p_x v_x p_y v_y;
x = [p_x; v_x];

A = zeros(2);
A(1, 2) = 1;

B = [0; 1];

f = A * x;
g = B;

x1min = -10;
x1mx = 10;
x2min = -10;
x2max = 10;

%% CLF and ECBF+
symbolic_clf = defineClf1(x,x1d);
symbolic_ecbf = defineECBF(x);

[x1,xdim,udim,f1,g1,clf,lf_clf,lg_clf,ecbf,dcbf1,lf_cbf1_1,lf_cbf1_2,lglf_cbf1_2,dcbf2,lf_cbf2_1,lf_cbf2_2,lglf_cbf2_2,lf_cbf3,lg_cbf3,lf_cbf4,lg_cbf4] = initSys(x, f, g, symbolic_ecbf, symbolic_clf);

%% 

odeFun = @dynamics; %system dynamics
controller = @ctrlCbfClfQp; %Exponential CBF and CLF controller


total_k = ceil(sim_t / dt); %total simulation time
x = x0; 
x_hold = x0;
t = 0; 

% % initialize traces.
xs = zeros(total_k, xdim);
xhold_s = zeros(total_k, xdim);
ts = zeros(total_k, 1);
ts_hold = zeros(total_k, 1);
Tclfs = []; 
Tcbfs = [];
tk_plusOne = [];
us = zeros(total_k-1, udim);
uhold_s = zeros(total_k-1, udim);
slacks = zeros(total_k-1, 1);
hs = zeros(total_k-1, 1);
Vs = zeros(total_k-1, 1);
xs(1, :) = x0';
xhold_s(1,:) = x0';
ts(1) = t;
% u_prev = [0;0];
Ts = 0.1; % sampling period (adjust according to your requirements)
% counter = 0;


current_timestep = [];
fk = [];
gk = [];
po = [];
u_hold = 0;
tnext = 0;
for k = 1:total_k-1
    if t + dt >= tnext
      [u, slack, h, V] = controller(x,clf,lf_clf,lg_clf,ecbf,dcbf1,lf_cbf1_1,lf_cbf1_2,lglf_cbf1_2,dcbf2,lf_cbf2_1,lf_cbf2_2,lglf_cbf2_2,lf_cbf3,lg_cbf3,lf_cbf4,lg_cbf4,udim);    
        us_temp = u;
        tk = ts(k);  % This is the current time step
        [roots_vector_secant, roots_vector_Bisection, Tcbf] = CBF_safe_period(x, ecbf, f1, g1, u, tk);
        all_roots_vector_secant = [all_roots_vector_secant roots_vector_secant'];
        all_roots_vector_Bisection = [all_roots_vector_Bisection roots_vector_Bisection'];
    
        Tcbfs = [Tcbfs Tcbf];
        Tclf = CLF_safe_period(x, clf, u,x1d,x2d);
        Tclfs = [Tclfs Tclf];
        po = [po min(Tcbf,Tclf)];
%         tnext = tk + min(Tcbf,Tclf);
        tnext = tk + 0.5;
        t_next_temp = tnext;
        us(k, :) = u';
        tk_plusOne = [tk_plusOne, tnext];
    else
        u_hold  = us_temp;
        tnext = t_next_temp;
        tk_plusOne = [tk_plusOne, tnext];
        us(k, :) = u_hold';
    end
   [ts_temp, xs_temp] = ode45(@(t, s) odeFun(t, s, us(k,:),f1,g1), [t t+dt], x);
    x = xs_temp(end, :)';

    xs(k+1, :) = x';
   
    ts(k+1) = t +dt;
    t = t + dt;
end



%% Plotting
% 
figure(1)
grid on
grid minor
hold on
plot(ts(1:end-1), us)
xlabel('t')
ylabel('u')
title('Control u over time')
hold off
% 
% figure(2)
% grid on
% grid minor
% hold on
% plot(ts(1:end-1), Tclfs')
% xlabel('t')
% ylabel('tclf')
% title('tclf input vs time')
% hold off
% 
% figure(3)
% grid on
% grid minor
% hold on
% plot(ts(1:end-1), Tcbfs')
% xlabel('t')
% ylabel('tcbf')
% title('tcbf input vs time')
% hold off

% figure(4)
% grid on
% grid minor
% hold on
% plot(ts(1:end-1), us)
% xlabel('t')
% ylabel('us')
% title('u input vs time')
% hold off
% 
figure(5)
grid on
grid minor
hold on
plot(xs(:,1), xs(:,2))
axis([-10 10 -10 10])
xlabel('x2')
ylabel('x1')
title('Periodic controller trajectory with ZOH')
hold off
% 
% figure(6)
% grid on
% grid minor
% hold on
% plot(xhold_s(:,1), xhold_s(:,2))
% xlabel('t')
% ylabel('xhold')
% title('trajectory (ZOH)')
% hold off
% 
figure(7)
grid on
grid minor
hold on
plot(ts(1:end-1), tk_plusOne' - ts(1:end-1))
xlabel('t')
ylabel('tnext')
title('tnext  input vs time')
hold off

figure(8)
grid on
grid minor
hold on
plot(ts, xs(:,1))

xlabel('t')
ylabel('position')
title('position x1 over time')
hold off
figure(9)
grid on
grid minor
hold on
plot(ts, xs(:,2))

xlabel('t')
ylabel('velocity')
title('velocity x2 over time')
hold off


% end

%% Functions


function Tclf = CLF_safe_period(symbolic_state, clf, uk, x1d, x2d)

            x = symbolic_state;
            
%             x1d = -10;
%             x2d = 0;

            V = clf(x);

%             D = 2*V + 2 * abs(uk) * sqrt(V) + 3 * abs(sqrt(V)) * abs(uk) + 2*(abs(uk))^2;
%             Tclf = (-2*(2*(x(2)) * (x(1) - x1d) + x(2)^2 + ((x(1) - x1d) + 2*x(2))*uk)) / D;

            D = 2*V + 2 * abs(uk) * sqrt(V) + 3 * sqrt(V) * abs(uk) + 2*(abs(uk))^2;
            Tclf = (-2*(2*x(2)*(x(1) - x1d) + x(2)^2 + ((x(1) - x1d) + 2*x(2)) * uk)) / D;
            
end


function [roots_vector_secant, roots_vector_Bisection, Tcbf] = CBF_safe_period(symbolic_state, ecbf, f1, g1, u, tk)

    % Define x, rho, k1, k2, B, and uk as given in your question
    x = symbolic_state;
    B = ecbf(x);
    k1 = 105;
    k2 = 20.5;
    L = 1;

   

%     rho = ( ((norm(f1(x) + g1(x) * uk)) / L) * exp(L * (t - tk)) - ((1 / L) * (norm(f1(x) + g1(x) * uk))) ); % Define rho value
    uk = u; % Define uk value

    ro = 1/L * (norm(f1(x) + g1(x) * uk));
%     rho = ro * exp(L * (t - tk)) - ( 1/L * (norm(f1(x) + g1(x) * uk)) );
    
    % Define zetabar functions
    zetabar_1 = @(t) ( k1*(x(2) -  ( ro * exp(L * (t - tk)) - ( 1/L * (norm(f1(x) + g1(x) * uk)) ) )    ) - k2*abs(uk) ) * t + B(1,1);
    zetabar_2 = @(t) ( -k1*(x(2) + ( ro * exp(L * (t - tk)) - ( 1/L * (norm(f1(x) + g1(x) * uk)) ) )    ) - k2*abs(uk) ) * t + B(2,1);
    zetabar_3 = @(t) -k1 * abs(uk) * t + B(3,1);
    zetabar_4 = @(t) -k1 * abs(uk) * t + B(4,1);
    
    % Define the secant method function

    % Find the roots for zetabar equations using the secant method
    tol = 0.001; % Tolerance for the secant method
    max_iter = 100; % Maximum number of iterations
    %%%%%%%%%%%%%%%initial guesses%%%%%%%%%%%%%%%%%%%%%%%%
    x0_1 = 0;
    x1_1 = 10;

    x0_2 = 0;
    x1_2 = 10;

    x0_3 = 0;
    x1_3 = 10;

    x0_4 = 0; % Initial guess 1
    x1_4 = 10;  % Initial guess 2
    %%%%%%%%%%%%%%%initial guesses%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%secant%%%%%
%     root_zetabar_1 = secant_method(zetabar_1, x0_1, x1_1, max_iter);
%     root_zetabar_2 = secant_method(zetabar_2, x0_2, x1_2, max_iter);
%     root_zetabar_3 = secant_method(zetabar_3, x0_3, x1_3, max_iter);
%     root_zetabar_4 = secant_method(zetabar_4, x0_4, x1_4, max_iter);
% 
%     roots_vector_secant = [root_zetabar_1,root_zetabar_2,root_zetabar_3,root_zetabar_4];
    roots_vector_secant = [0,0,0,0];

    %%%%Bisection%%%%
    a1 = 0;
    b1 = 40;

    a2 = 0;
    b2 = 45;

    a3 = 0;
    b3 = 1e+12;

    a4 = 0;
    b4 = 1e+12;

    root_zetabar_B1 = bisection(zetabar_1, a1, b1, tol);
    root_zetabar_B2 = bisection(zetabar_2, a2, b2, tol);
    root_zetabar_B3 = bisection(zetabar_3, a3, b3, tol);
    root_zetabar_B4 = bisection(zetabar_4, a4, b4, tol);

    roots_vector_Bisection = [root_zetabar_B1, root_zetabar_B2, root_zetabar_B3, root_zetabar_B4];
%     roots_vector_Bisection = [0,0,0,0];

    %%%%%%%%%%%%%%%%%
    
    Tcbf = min(roots_vector_Bisection);

end

%%%%%%%%%%%%%%%%%%%%%%%Finding roots fucntions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function root = secant_method(func, x0, x1, tol, max_iter)
%     for i = 1:max_iter
%         x2 = x1 - func(x1) * (x1 - x0) / (func(x1) - func(x0));
%         if abs(x2 - x1) < tol
%             root = x2;
%             return
%         end
%         x0 = x1;
%         x1 = x2;
%     end
%     root = x1;
%     fprintf('Maximum iterations reached. The current approximation is %f\n', root);
% end

function y = secant_method(f,a,b,maxerr)

c = (a*f(b) - b*f(a))/(f(b) - f(a));
disp('   Xn-1      f(Xn-1)      Xn      f(Xn)      Xn+1      f(Xn+1)');
disp([a f(a) b f(b) c f(c)]);
while abs(f(c)) > maxerr
    a = b;
    b = c;
    c = (a*f(b) - b*f(a))/(f(b) - f(a));
    disp([a f(a) b f(b) c f(c)]);

    flag = flag + 1;

    if(flag == 100)
        break;
    end
end

display(['Root is x = ' num2str(c)]);
y = c;
end

% function root = Bisection_method(func, a, b, tol)
%     while (b - a) > tol
%     c = (a + b) / 2;
%         if func(a) * func(c) < 0 %checks if root of function lies within the interval [a,c]
%             b = c;
%         else
%             a = c;
%         end
%     end
%     root = (a + b) / 2;
% end

function  m = bisection(f, low, high, tol)
disp('Bisection Method'); 

% Evaluate both ends of the interval
y1 = feval(f, low);
y2 = feval(f, high);
i = 0; 

% Display error and finish if signs are not different
if y1 * y2 > 0
   disp('Have not found a change in sign. Will not continue...');
   m = 6666666
   return
end 

% Work with the limits modifying them until you find
% a function close enough to zero.
disp('Iter    low        high          x0');
    while (abs(high - low) >= tol)
        i = i + 1;
        % Find a new value to be tested as a root
        m = (high + low)/2;
        y3 = feval(f, m);
        if y3 == 0
            fprintf('Root at x = %f \n\n', m);
            return
        end
        fprintf('%2i \t %f \t %f \t %f \n', i-1, low, high, m);   
    
        % Update the limits
        if y1 * y3 > 0
            low = m;
            y1 = y3;
        else
            high = m;
        end
    end 

end

% Show 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dx = dynamics(t, x, u,f1,g1)

            dx = f1(x) + g1(x) * u;
end


function clf = defineClf1(symbolic_state,x1d)

            x = symbolic_state;
            clf = [x(1) - x1d; x(2)]' * [1, 0.5 ; 0.5 ,1] * [ x(1) - x1d; x(2)];        
end


function ecbf = defineECBF(symbolic_state)
            
            x1min = -10;
            x1max = 10;
            x2min = -10;
            x2max = 10;
            x = symbolic_state;

            h1 = x(1) - x1min;
            h2 = -x(1) + x1max;
            h3 = x(2) - x2min;
            h4 = -x(2) + x2max;

            ecbf = [h1;h2;h3;h4];
    
end
                

function [x1,xdim,udim,f1,g1,clf,lf_clf,lg_clf,ecbf,dcbf1,lf_cbf1_1,lf_cbf1_2,lglf_cbf1_2,dcbf2,lf_cbf2_1,lf_cbf2_2,lglf_cbf2_2,lf_cbf3,lg_cbf3,lf_cbf4,lg_cbf4] = initSys(symbolic_x, symbolic_f, symbolic_g, symbolic_ecbf, symbolic_clf)
    if isempty(symbolic_x) || isempty(symbolic_f) || isempty(symbolic_g)
        error('x, f, g is empty.');
    end

    if ~isa(symbolic_f, 'sym')
        f_ = sym(symbolic_f);
    else
        f_ = symbolic_f;
    end
    if ~isa(symbolic_g, 'sym')
        g_ = sym(symbolic_g);
    else
        g_ = symbolic_g;
    end

    x1 = symbolic_x;
    % Setting state and input dimension.
    xdim = size(x1, 1);
    udim = size(g_, 2);

    % Setting f and g (dynamics)
    f1 = matlabFunction(f_, 'vars', {x1});
    g1 = matlabFunction(g_, 'vars', {x1});            

    % Obtaining Lie derivatives of CBF.
    if ~isempty(symbolic_ecbf)
        fprintf('n')
        
        %h1(x)
        dcbf1_ = simplify(jacobian(symbolic_ecbf(1,1), x1));
        lf_cbf1_1_ = dcbf1_ * f_;
        ddcbf1 = simplify(jacobian(lf_cbf1_1_, x1));
        lf_cbf1_2_ = ddcbf1 * f_;
        lglf_cbf1_2_ = ddcbf1 * g_;
        
        %h2(x)
        dcbf2_ = simplify(jacobian(symbolic_ecbf(2,1), x1));
        lf_cbf2_1_ = dcbf2_ * f_;
        ddcbf2 = simplify(jacobian(lf_cbf2_1_, x1));
        lf_cbf2_2_ = ddcbf2 * f_;
        lglf_cbf2_2_ = ddcbf2 * g_;

        %h3(x)
        dcbf3 = simplify(jacobian(symbolic_ecbf(3,1), x1));
        lf_cbf3_ = dcbf3 * f_;
        lg_cbf3_ = dcbf3 * g_;
        
        %h4(x)
        dcbf4 = simplify(jacobian(symbolic_ecbf(4,1), x1));
        lf_cbf4_ = dcbf4 * f_;
        lg_cbf4_ = dcbf4 * g_;

        %%%%%%%%%%%%%%%%%matlabfunction%%%%%%%%%%%%%%%%%%%%%
        ecbf = matlabFunction(symbolic_ecbf, 'vars', {x1});

        %h1(x)
        dcbf1 = matlabFunction(dcbf1_, 'vars', {x1});
        lf_cbf1_1 = matlabFunction(lf_cbf1_1_, 'vars', {x1});
        lf_cbf1_2 = matlabFunction(lf_cbf1_2_, 'vars', {x1});
        lglf_cbf1_2 = matlabFunction(lglf_cbf1_2_, 'vars', {x1});

        %h2(x)
        dcbf2 = matlabFunction(dcbf2_, 'vars', {x1});
        lf_cbf2_1 = matlabFunction(lf_cbf2_1_, 'vars', {x1});
        lf_cbf2_2 = matlabFunction(lf_cbf2_2_, 'vars', {x1});
        lglf_cbf2_2 = matlabFunction(lglf_cbf2_2_, 'vars', {x1});

        %h3(x)
        lf_cbf3 = matlabFunction(lf_cbf3_, 'vars', {x1});
        lg_cbf3 = matlabFunction(lg_cbf3_, 'vars', {x1});

        %h4(x)
        lf_cbf4 = matlabFunction(lf_cbf4_, 'vars', {x1});
        lg_cbf4 = matlabFunction(lg_cbf4_, 'vars', {x1});
     
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    % Obtaining Lie derivatives of CLF.    
    if ~isempty(symbolic_clf)
        dclf = simplify(jacobian(symbolic_clf, x1)); %computes deltaV
        lf_clf_ = dclf * f_;
        lg_clf_ = dclf * g_;
        clf = matlabFunction(symbolic_clf, 'vars', {x1});                       
        lf_clf = matlabFunction(lf_clf_, 'vars', {x1});
        lg_clf = matlabFunction(lg_clf_, 'vars', {x1});        
    end
end

%% controller
function [u, slack, B, V, feas, comp_time] = ctrlCbfClfQp(x,clf,lf_clf,lg_clf,ecbf,dcbf1,lf_cbf1_1,lf_cbf1_2,lglf_cbf1_2,dcbf2,lf_cbf2_1,lf_cbf2_2,lglf_cbf2_2,lf_cbf3,lg_cbf3,lf_cbf4,lg_cbf4,udim)
    

    clf_rate = 0.8; %epsilon
    cbf_rate = [105, 20.5]; %Kb
    
    u_max = 20;
    u_min  = -20;
    
    weight_slack = 1;
    weight_input = 5;
               
    tstart = tic;
    V = clf(x);
    LfV = lf_clf(x);
    LgV = lg_clf(x);

    B = ecbf(x);
   
    %h1(x)
    h1_dotx = dcbf1(x);
    Lf1_h1 = lf_cbf1_1(x);
    Lf2_h1 = lf_cbf1_2(x);
    LgLf_h1 = lglf_cbf1_2(x);


    %h2(x)
    h2_dotx = dcbf2(x);
    Lf1_h2 = lf_cbf2_1(x);
    Lf2_h2 = lf_cbf2_2(x);
    LgLf_h2 = lglf_cbf2_2(x);


    %h3(x)
    Lf_h3 = lf_cbf3(x);
    Lg_h3 = lg_cbf3(x);

    %h4(x)
    Lf_h4 = lf_cbf4(x);
    Lg_h4 = lg_cbf4(x);


    u_ref = zeros(udim, 1);
    with_slack = 1;


     if size(u_ref, 1) ~= udim
        error("Wrong size of u_ref, it should be (udim, 1) array.");
    end                
        
    %% Constraints: A[u; slack] <= b adds inp
    if with_slack
        % CLF and CBF constraints.
        A = [LgV, -1;
            -LgLf_h1, 0;
            -LgLf_h2, 0;
            -Lg_h3, 0;
            -Lg_h4, 0];
        A(1,2) = 0;
        b = [-LfV - clf_rate * V;
            Lf2_h1 - cbf_rate(1,1)*B(1,1)*Lf1_h1 - cbf_rate(1,2)*(Lf1_h1 - cbf_rate(1,1)*B(1,1));
            Lf1_h2 - cbf_rate(1,1)*B(2,1)*Lf1_h2 - cbf_rate(1,2)*(Lf1_h2 - cbf_rate(1,1)*B(2,1));
            Lf_h3 + cbf_rate(1,1) * B(3,1);
            Lf_h4 + cbf_rate(1,1) * B(4,1)]; 

        szA1 = size(A);
        szB1 = size(b);
        % Add input constraints if u_max or u_min exists.
        if exist('u_max', 'var')
            A = [A; eye(udim), zeros(udim, 1);];
            b = [b; u_max * ones(udim, 1)];
        end
        if exist('u_min', 'var')
            A = [A; -eye(udim), zeros(udim, 1);];
            b = [b; -u_min * ones(udim, 1)];
        end        

    end

    %% Cost
    if exist('weight_input', 'var')
            weight_input1 = weight_input * eye(udim);

    else
        weight_input1 = eye(udim);
    end

    if with_slack     

        H = [weight_input1, zeros(udim, 1);
            zeros(1, udim), weight_slack];
        f_ = [-weight_input1 * u_ref; 0];
        [u_slack, ~, exitflag, ~] = quadprog(H, f_, A, b, [], [], [], [], [],optimset('Display', 'off'));

        
        if exitflag == -2 
            feas = 0;
            disp("Infeasible QP. CBF constraint is conflicting with input constraints.");
        else
            feas = 1;
        end
        u = u_slack(1:udim);
        slack = u_slack(end);

    end
    comp_time = toc(tstart);
end