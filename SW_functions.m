function f_ts = f_ftcs(ts,f_0,idim,del_x,del_t,U_mean) 
%% F_FTCS Function
% ftcs, height distribution after ts time steps
% forward Euler time extrapolation & center difference
% ts :  total time steps
% h_0 : height distribution at initial state
% h_ts : height distribution after ts time steps
f_ts = zeros(ts+1,idim); % height distribution array, each row is h(x) at time t
f_ts(1,:) = f_0; % t=0 at initial state.
for t = 1:ts
    dfdt = fTend(f_ts(t,:),del_x,U_mean); % fTend is defined later, center differnce
    f_ts(t+1,:) = t_forewardE(f_ts(t,:),del_t,dfdt); % t_LeapFrog is defined later
    f_ts(t+1,:) = BC_periodic(f_ts(t+1,:),idim); % BC_periodic is defined later
end 
end
% ctcs, height distribution after ts time steps
% leap-frog time extraolation & center difference
function f_ts = f_ctcs(ts,f_0,idim,del_x,del_t,U_mean)
% ts :  total time steps
% h_0 : height distribution at initial state
% h_ts : height distribution after ts time steps
f_ts = zeros(ts+1,idim); % height distribution array, each row is h(x) at time t
f_ts(1,:) = f_0; % t=0 at initial state.
% step 1
dfdt = fTend(f_ts(1,:),del_x,U_mean); % fTend is defined later
f_ts(2,:) = t_forewardE(f_ts(1,:),del_t,dfdt); % t_forewardE is defined later
f_ts(2,:) = BC_periodic(f_ts(2,:),idim); % BC_periodic is defined later
% step 2 and later, until step ts, at ts+1 row
for t = 2:ts
    dfdt = fTend(f_ts(t,:),del_x,U_mean); % fTend is defined later
    f_ts(t+1,:) = t_LeapFrog(f_ts(t-1,:),del_t,dfdt); % t_LeapFrog is defined later
    f_ts(t+1,:) = BC_periodic(f_ts(t+1,:),idim); % BC_periodic is defined later
end 
end
% periodic boundary condition
function f = BC_periodic(f,idim)
% f input : without boundary elements
% f output : all elements
    f(1) = f(idim-1);
    f(idim) = f(2);
end
% forward Euler method, time extrapolation
function fp = t_forewardE(f,del_t,fTend) 
% fp : f at t+1
% f : f at t
% fTend : f tendency at t
% del_T : time step, s
    
    fp = f + del_t*fTend;
end
% leap-frog method, time extrapolation
function fp = t_LeapFrog(fm,del_t,fTend) 
% fp : f at t+1
% fm : f at t-1
% fTend : f tendency at t
% del_T : time step, s
    
    fp = fm + 2*del_t*fTend;
end
% tendency
function dfdt = fTend(f,del_x,U_mean) 
% dfdt :  f tendency
% f : input f distribution at anytime
dfdx = gradient(f,del_x); 
dfdx(1) = 0;
dfdx(end) = 0;
dfdt = -U_mean*dfdx; % ∂h/∂t = - mean(U) ∂h/∂x
end
% Hovmoeller diagram
function hov = hovgram(ft,x,y)
hov = contourf(x,y,ft,'linecolor','none');
title("Hovmoeller diagram of height");
xlabel("distance-x (m)");
ylabel("time-t (s)");
c = colorbar;
c.Label.String = 'height (m)';
end