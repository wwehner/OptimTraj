function [t,x,u,defects,pathCost,dxdalpha,dJdalpha] = multiShootIntegrate(intScheme, tSpan, state, control, pack, dynFun, pathObj, gradInfo)
% this function carries out the integration for the Multiple Shooting
% approach. Several different integration approaches are available defined
% by intScheme
% intScheme = 1 (classical rungekutta, explicit)
%           = 2 (backward euler, implicit)
%           = 3 (trapezoidal, implicit)
%           = 4 (hermite-Simpson, implicit)

nState = pack.nState;
nControl = pack.nControl;
nSegment = pack.nSegment;
nSubStep = pack.nSubStep;

% options for fsolve (used for impicit integration schemes)
options = optimset('Display','off','Jacobian','on');

% NOTES:
%   The following bit of code is a bit confusing, mostly due to the
%   need for vectorization to make things run at a reasonable speed in
%   Matlab. Part of the confusion comes because the decision variables
%   include the state at the beginning of each segment, but the control
%   at the beginning and middle of each substep - thus there are more
%   control grid-points than state grid points. The calculations are
%   vectorized over segments, but not sub-steps, since the result of
%   one sub-step is required for the next.

% time, state, and control at the ends of each substep
nTime = 1+nSegment*nSubStep;
t = linspace(tSpan(1), tSpan(2), nTime);
x = zeros(nState, nTime);
u = control(:,1:2:end); % Control a the endpoints of each segment
uMid = control(:,2:2:end);  %Control at the mid-points of each segment
c = zeros(1, nTime-1);  %Integral cost for each segment
dt = (t(end)-t(1))/(nTime-1);

idx = 1:nSubStep:(nTime-1);   %Indicies for the start of each segment
x(:,[idx,end]) = state;   %Fill in the states that we already know

% VARIABLES for analytic gradient evaluations.
% let alpha = decVars
% alpha = [t0, tF, x0, x1, ..., xN, u0, uM0, u1, ..., uN]
% size alpha = 2 (t0, tF) + nstate*(nSegment+1) + nControl* ...
% dxdalpha = partial derivative of state w.r.t. alpha
nalpha = 2 + nState*(1+nSegment) + nControl*(1+2*nSubStep*nSegment);

% need to calculate dxdalpha separately for each segment of the
% trajectory
dxdalpha = cell(1,nSegment);
for i = 1:nSegment
    dxdalpha{i} = zeros(nState,nalpha,nSubStep+1);
    cols = 2+(i-1)*nState+(1:nState);

    % since the state at the start of each segment is included as a
    % alpha, dxdalpha[0] = I
    dxdalpha{i}(:,cols,1) = eye(nState);
end

% T = total time of trajectory, the gradient of T is included since t0,
% and tF are included as alpha
dTdalpha = zeros(1,nalpha); dTdalpha(1:2) = [-1,1];

% dt = time step, dt = (tF-t0)/N so it too has a partial derivative
% w.r.t. alpha
dt_dalpha = zeros(1,nalpha);

% t_n = discrete time variable
t_n = 0:nTime-1;

% gradient of path cost
dJdalpha = zeros(1,nalpha);

for iSubStep = 1:nSubStep
    % March forward Runge-Kutta step

    t0 = t(idx);
    x0 = x(:,idx);



    %------------------------------------------
    % Code for calculating dxdalpha (partial derivative of state w.r.t.
    % the descision parameters): dxdalpha = nstate x nalpha

    % Gradient of time w.r.t. decVars
    % ------------------------------------------------------------
    % dt = (tF-t0)/(nTime-1)
    % t = t0 + n*dt
    % t = t0 + n*(tF-t0)/(nTime-1)
    % t = t0*(1-n/(nTime-1)) + tF*(n/(nTime-1))
    %
    % alpha = [t0, tF, x0, x1, ..., xN, u0, uM0, u1, ..., uN]
    % dt/dalpha = [1 - n/(nTime-1), n/(nTime-1), 0, 0, ... 0]
    % ------------------------------------------------------------

    % Switch over possible integration schemes
    % rungeKutta = 1, trapezoidal = 2, hermite-simpson = 3
    switch intScheme

        % -------------------------------------------------------------
        % RungeKutta Integration

        case 1 % rungeKutta

            [k0, dk0] = combinedDynGrad(t0,        x0,                         u(:,idx), dynFun,pathObj);
            [k1, dk1] = combinedDynGrad(t0+0.5*dt, x0 + 0.5*dt*k0(1:nState,:), uMid(:,idx), dynFun,pathObj);
            [k2, dk2] = combinedDynGrad(t0+0.5*dt, x0 + 0.5*dt*k1(1:nState,:), uMid(:,idx), dynFun,pathObj);
            [k3, dk3] = combinedDynGrad(t0+dt,     x0 +     dt*k2(1:nState,:), u(:,idx+1), dynFun,pathObj);
            z = (dt/6)*(k0 + 2*k1 + 2*k2 + k3);  %Change over the sub-step
            
            % discrete time on each segment
            t_n0 = t_n(idx);

            for j = 1:nSegment

                % d(t[n])/dalpha
                dt_dalpha(1) = 1 - t_n0(j)/(nTime-1);
                dt_dalpha(2) = t_n0(j)/(nTime-1);

                % du[n]/dalpha
                du_dalpha = zeros(nControl,nalpha);
                du_dalpha(:,gradInfo.indu(:,idx(j))) = eye(nControl);

                % duMid[n]/dalpha
                duMid_dalpha = zeros(nControl,nalpha);
                duMid_dalpha(:,gradInfo.indumid(:,idx(j))) = eye(nControl);

                % du[n+1]/dalpha
                du1_dalpha = zeros(nControl,nalpha);
                du1_dalpha(:,gradInfo.indu(:,idx(j)+1)) = eye(nControl);

                % dk0/dalpha
                dk0da = dk0(:,:,j) * [dt_dalpha; dxdalpha{j}(:,:,iSubStep); du_dalpha];

                % dk1/dalpha
                dk1da = dk1(:,:,j) * [dt_dalpha + 0.5/(nTime-1)*dTdalpha; dxdalpha{j}(:,:,iSubStep) + 0.5*dt*dk0da(1:nState,:) + 0.5/(nTime-1)*k0(1:nState,j)*dTdalpha; duMid_dalpha];

                % dk2/dalpha
                dk2da = dk2(:,:,j) * [dt_dalpha + 0.5/(nTime-1)*dTdalpha; dxdalpha{j}(:,:,iSubStep) + 0.5*dt*dk1da(1:nState,:) + 0.5/(nTime-1)*k1(1:nState,j)*dTdalpha; duMid_dalpha];

                % dk3/dalpha
                dk3da = dk3(:,:,j) * [dt_dalpha + 1/(nTime-1)*dTdalpha; dxdalpha{j}(:,:,iSubStep) + dt*dk2da(1:nState,:) + 1/(nTime-1)*k2(1:nState,j)*dTdalpha; du1_dalpha];

                dz = (dt/6)*(dk0da + 2*dk1da + 2*dk2da + dk3da)...
                    + 1/(6*(nTime-1))*(k0(:,j)+2*k1(:,j)+2*k2(:,j)+k3(:,j))*dTdalpha;

                % update dxdalpha
                dxdalpha{j}(:,:,iSubStep+1) = dxdalpha{j}(:,:,iSubStep) + dz(1:nState,:);

                % update dJdalpha
                dJdalpha  = dJdalpha + dz(nState+1,:);
            end

        % -------------------------------------------------------------
        % Backward Euler Integration
        %     x[n+1] = x[n] + dt * f(x[n+1],u[n+1])
        
        case 2 

            % z = x[n+1] - x[n]
            z = zeros(nState+1,nSegment);
            
            
            % % for loop integration with fsolve, slow but easier to understand
            %{
            for i_test = 1:nSegment
                z(:,i_test) = fsolve(@(x) backwardEuler(x, x0(:,i_test),...
                    u(:,idx(i_test)), u(:,idx(i_test)+1), dynFun, pathObj,...
                    t0(i_test),dt), z(:,i_test), options);
            end
            %}
            
            % vectorized fsolve (same as block above)
            z = fsolve(@(x) backwardEuler(x, x0, u(:,idx), u(:,idx+1), ...
                dynFun, pathObj, t0,dt), z(:), options);
            z = reshape(z,nState+1,length(t0));

            % state and control at next time step (xn1 = x[n+1], un1 = u[n+1])
            xn1 = x0 + z(1:nState,:);
            un1 = u(:,idx+1);

            % fn1 = f(x[n+1],u[n+1]), dfn1 = jacobian w.r.t. (t,x,u)
            [fn1,dfn1] = combinedDynGrad(t0+dt, xn1, un1, dynFun, pathObj);
            
            % discrete time on each segment
            t_n0 = t_n(idx);

            % Calculate dx/dalpha (gradient of state w.r.t. alpha)
            for j = 1:nSegment
                
                M_alpha = eye(nState);
                
                % d(t[n])/dalpha
                dt_dalpha(1) = 1 - t_n0(j)/(nTime-1);
                dt_dalpha(2) = t_n0(j)/(nTime-1);

                % d(u[n+1])/dalpha
                dudalpha = zeros(nControl,nalpha);
                cols = gradInfo.indu(:,idx(j)+1);
                dudalpha(:,cols) = eye(nControl);

                % dfn1/dalpha (only include part on control)
                dfn1du_duda = dfn1(1:end-1,[1 2+nState:end],j) * [dt_dalpha + 1/(nTime-1)*dTdalpha; dudalpha];

                M_alpha = M_alpha - dt*(dfn1(1:end-1,2:nState+1,j));

                dxdalpha{j}(:,:,iSubStep+1) = (M_alpha) \ (dxdalpha{j}(:,:,iSubStep) + (dt)*(dfn1du_duda)...
                + 1/((nTime-1))*(fn1(1:end-1,j))*dTdalpha);

                dJdalpha = dJdalpha + 1/(nTime-1)*fn1(end,j)*dTdalpha + ...
                    dt*(dfn1(end,:,j) * [dt_dalpha + 1/(nTime-1)*dTdalpha; dxdalpha{j}(:,:,iSubStep+1); dudalpha]);

            end

        % -------------------------------------------------------------
        % Trapezoidal Integration
        %     x[n+1] = x[n] + dt/2 * ( f(x[n+1],u[n+1]) + f(x[n],u[n]) )
        
        case 3
            
            % z = x[n+1] - x[n]
            z = zeros(nState+1,nSegment);

            % vectorized fsolve (same as block above)
            z = fsolve(@(x) trapezoidal(x,x0,u(:,idx),u(:,idx+1),...
                dynFun,pathObj,t0,dt),z(:),options);
            z = reshape(z,nState+1,length(t0));
            
            % state and control at next time step (xn1 = x[n+1], un1 = u[n+1])
            xn = x0;
            xn1 = x0+z(1:nState,:);
            un = u(:,idx);
            un1 = u(:,idx+1);

            % fn1 = f(x[n+1],u[n+1]), dfn1 = jacobian w.r.t. (t,x,u)
            [fk,dfk] = combinedDynGrad(t0, xn, un, dynFun, pathObj);
            [fk1,dfk1] = combinedDynGrad(t0+dt, xn1, un1, dynFun, pathObj);
            
            % discrete time on each segment
            t_n0 = t_n(idx);

            % Calculate dx/dalpha (gradient of state w.r.t. alpha)
            for j = 1:nSegment
                
                M_alpha = eye(nState);
                
                % d(t[n])/dalpha
                dt_dalpha(1) = 1 - t_n0(j)/(nTime-1);
                dt_dalpha(2) = t_n0(j)/(nTime-1);

                % dfk/dalpha
                dudalpha = zeros(nControl,nalpha);
                cols = gradInfo.indu(:,idx(j));
                dudalpha(:,cols) = eye(nControl);
                dfkda = dfk(:,:,j) * [dt_dalpha; dxdalpha{j}(:,:,iSubStep); dudalpha];

                % dfk1/dalpha (only include part on control)
                dudalpha = zeros(nControl,nalpha);
                cols = gradInfo.indu(:,idx(j)+1);
                dudalpha(:,cols) = eye(nControl);
                dfk1du_duda = dfk1(:,[1 2+nState:end],j) * [dt_dalpha + 1/(nTime-1)*dTdalpha; dudalpha];

                M_alpha = M_alpha - dt/2*(dfk1(1:end-1,2:nState+1,j));

                dxdalpha{j}(:,:,iSubStep+1) = (M_alpha) \ (dxdalpha{j}(:,:,iSubStep) + (dt/2)*(dfkda(1:end-1,:) + dfk1du_duda(1:end-1,:))...
                  + 1/(2*(nTime-1))*(fk(1:end-1,j)+fk1(1:end-1,j))*dTdalpha);
              
                % dJ_dalpha
                dJdalpha = dJdalpha + 1/(2*(nTime-1))*(fk(end,j)+fk1(end,j))*dTdalpha + ...
                    dt/2 * (dfkda(end,:) + dfk1(end,:,j)*[dt_dalpha + 1/(nTime-1)*dTdalpha; dxdalpha{j}(:,:,iSubStep+1);dudalpha]);

            end

        % -------------------------------------------------------------
        % Hermite-simpson Integration

        case 4 

            % z = x[n+1] - x[n]
            z = zeros(nState+1,nSegment);
            
            % for loop integration with fsolve, slow but easier to understand
            %{
            for i_test = 1:nSegment
              z(:,i_test) = fsolve(@(x) hermiteSimpson(x,x0(:,i_test),...
                  u(:,idx(i_test)),uMid(:,idx(i_test)),u(:,idx(i_test)+1),...
                  dynFun,pathObj,t0(i_test),dt),z(:,i_test),options);
            end
            %}
            
            % vectorized fsolve (same as block above)
            z = fsolve(@(x) hermiteSimpson_v(x,x0,u(:,idx),uMid(:,idx),u(:,idx+1),...
              dynFun,pathObj,t0,dt),z(:),options);
            z = reshape(z,nState+1,length(t0));

            % state and control at next time step (xn1 = x[n+1], un1 = u[n+1])
            xn = x0;
            xn1 = x0+z(1:nState,:);
            un = u(:,idx);
            un1 = u(:,idx+1);
            uM = uMid(:,idx);


            [fk,dfk] = combinedDynGrad(t0, xn, un, dynFun,pathObj);
            [fk1,dfk1] = combinedDynGrad(t0+dt, xn1, un1 ,dynFun,pathObj);

            ybar = 0.5 * (xn1 + xn) + dt/8 * (fk(1:end-1,:) - fk1(1:end-1,:));

            [fbar,dfbar] = combinedDynGrad(t0+0.5*dt, ybar, uM, dynFun,pathObj);

            % discrete time on each segment
            t_n0 = t_n(idx);
            
            % Calculate dx/dalpha (gradient of state w.r.t. alpha)
            for j = 1:nSegment
                
                M_alpha = eye(nState);

                % d(t[n])/dalpha
                dt_dalpha(1) = 1 - t_n0(j)/(nTime-1);
                dt_dalpha(2) = t_n0(j)/(nTime-1);

                % dfk/dalpha
                dudalpha = zeros(nControl,nalpha);
                cols = gradInfo.indu(:,idx(j));
                dudalpha(:,cols) = eye(nControl);
                dfkda = dfk(:,:,j) * [dt_dalpha; dxdalpha{j}(:,:,iSubStep); dudalpha];

                % dfk1/dalpha (only include part on control)
                dudalpha = zeros(nControl,nalpha);
                cols = gradInfo.indu(:,idx(j)+1);
                dudalpha(:,cols) = eye(nControl);
                dfk1du_duda = dfk1(:,[1 2+nState:end],j) * [dt_dalpha + 1/(nTime-1)*dTdalpha; dudalpha];

                M_alpha = M_alpha - dt/6*dfk1(1:end-1,2:nState+1,j);

                % dybar/dalpha
                dybarda = 0.5*dxdalpha{j}(:,:,iSubStep) + dt/8*(dfkda(1:end-1,:) - dfk1du_duda(1:end-1,:)) + ...
                1/(8*(nTime-1))*(fk(1:end-1,j)-fk1(1:end-1,j))*dTdalpha;

                % dfbar/dalpha
                dudalpha = zeros(nControl,nalpha);
                cols = gradInfo.indumid(:,idx(j));
                dudalpha(:,cols) = eye(nControl);
                dfbarda = dfbar(:,:,j) * [dt_dalpha + 0.5/(nTime-1)*dTdalpha; dybarda; dudalpha];

                M_alpha = M_alpha - dt/6*(4*dfbar(1:end-1,2:nState+1,j)*(.5*eye(nState) - dt/8*dfk1(1:end-1,2:nState+1,j)));


                dxdalpha{j}(:,:,iSubStep+1) = (M_alpha) \ (dxdalpha{j}(:,:,iSubStep) + (dt/6)*(dfkda(1:end-1,:) + 4*dfbarda(1:end-1,:) + dfk1du_duda(1:end-1,:))...
                + 1/(6*(nTime-1))*(fk(1:end-1,j)+4*fbar(1:end-1,j)+fk1(1:end-1,j))*dTdalpha);

                % x[k+1] = x[k] + dt/6 ( f[k] + 4fbar + f[k+1])

                % duM_dalpha
                duMdalpha = zeros(nControl,nalpha);
                cols = gradInfo.indumid(:,idx(j));
                duMdalpha(:,cols) = eye(nControl);

                % duk1_dalpha
                duk1dalpha = zeros(nControl,nalpha);
                cols = gradInfo.indu(:,idx(j)+1);
                duk1dalpha(:,cols) = eye(nControl);

                % dybar_dalpha
                dybarda = 0.5*(dxdalpha{j}(:,:,iSubStep+1) + dxdalpha{j}(:,:,iSubStep)) + ...
                dt/8*(dfkda(1:end-1,:)-dfk1(1:end-1,:,j)*[dt_dalpha + 1/(nTime-1)*dTdalpha; dxdalpha{j}(:,:,iSubStep+1); duk1dalpha]) +...
                1/(8*(nTime-1))*(fk(1:end-1,j)-fk1(1:end-1,j))*dTdalpha;


                % dJ_dalpha
                dJdalpha = dJdalpha + 1/(6*(nTime-1))*(fk(end,j)+4*fbar(end,j)+fk1(end,j))*dTdalpha + ...
                dt/6 * (dfkda(end,:) + 4*dfbar(end,:,j)*[dt_dalpha + 0.5/(nTime-1)*dTdalpha; dybarda; duMdalpha] + dfk1(end,:,j)*[dt_dalpha + 1/(nTime-1)*dTdalpha; dxdalpha{j}(:,:,iSubStep+1); duk1dalpha]);
            end

        otherwise
            error('Undefined Integration Method')

    end % end switch


    xNext = x0 + z(1:nState,:);  %Next state
    c(idx) = z(end,:);  %Integral of the cost function over this step

    if iSubStep == nSubStep %We've reached the end of the interval
        % Compute the defect vector:
        defects = xNext - x(:,idx+1);
    else
        % Store the state for next step in time
        idx = idx+1;   %  <-- This is important!!
        x(:,idx) = xNext;
    end


end

pathCost = sum(c);  %Sum up the integral cost over each segment


end

%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%%

function [dz, J] = combinedDynGrad(t,x,u,dynamics,pathObj)
% [dz, dJ] = combinedDynamics(t,x,u,dynamics,pathObj)
%
% This function packages the dynamics and the cost function together so
% that they can be integrated at the same time.
%
% INPUTS:
%   t = [1, nTime] = time vector (grid points)
%   x = [nState, nTime] = state vector at each grid point
%   u = [nControl, nTime] = control vector at each grid point
%   dynamics(t,x,u) = dynamics function handle
%               dx = [nState, nTime] = dx/dt = derivative of state wrt time
%   pathObj(t,x,u) = integral cost function handle
%                 dObj = [1, nTime] = integrand from the cost function
%
% OUTPUTS:
%   dz = [dx; dObj] = combined dynamics of state and cost
%   dJ = [JAC(dynamics), JAC(objective)] = combined jacobian of dynamics
%   and objective w.r.t. (t,x,u)

if nargout < 2
  
  dx = dynamics(t,x,u);
  if isempty(pathObj)
      dc = zeros(size(t));
  else
      dc = pathObj(t,x,u);
  end

  dz = [dx;dc];  %Combine and return

else
  
  nState = size(x,1);
  nControl = size(u,1);

  [dx,Jx] = dynamics(t,x,u);
  if isempty(pathObj)
      dc = zeros(size(t));
      Jc = zeros(1,1+nState+nControl,length(t));
  else
      [dc,Jc] = pathObj(t,x,u);
      Jc = reshape(Jc,1,1+nState+nControl,length(t));
  end

  dz = [dx;dc];

  J = cat(1,Jx,Jc);
end

end

function [dz, J] = combinedJacobian(t,x,u,dynamics,pathObj)
% [dz, dJ] = combinedDynamics(t,x,u,dynamics,pathObj)
%
% This function packages the dynamics and the cost function together so
% that they can be integrated at the same time.
%
% INPUTS:
%   t = [1, nTime] = time vector (grid points)
%   x = [nState, nTime] = state vector at each grid point
%   u = [nControl, nTime] = control vector at each grid point
%   dynamics(t,x,u) = dynamics function handle
%               dx = [nState, nTime] = dx/dt = derivative of state wrt time
%   pathObj(t,x,u) = integral cost function handle
%                 dObj = [1, nTime] = integrand from the cost function
%
% OUTPUTS:
%   dz = [dx; dObj] = combined dynamics of state and cost
%   dJ = [JAC(dynamics), JAC(objective)] = combined jacobian of dynamics
%   and objective w.r.t. (t,x,u)

if nargout < 2
  
    dx = dynamics(t,x,u);
    if isempty(pathObj)
      dc = zeros(size(t));
    else
      dc = pathObj(t,x,u);
    end

    dz = [dx;dc];  %Combine and return

else
  
    nTime = length(t);
    
    % Here, nStateCmbn = number of states + 1 to account for pathCost considered
    % as an additional state
    nState = size(x,1); 
    nStateCmbn = nState + 1;

    % Get dynamics and jacobian of dynamics which must be augmented with a
    % column of zeros to account for pathCost, cosidered as an additional
    % state
    [dx,Jx] = dynamics(t,x,u);

    % Jacobian w.r.t. x
    Jx = [Jx(:,1+(1:nState),:),zeros(nState,1,nTime)];

    if isempty(pathObj)
        dc = zeros(size(t));
        Jc = zeros(1,nStateCmbn,length(t));
    else
        [dc,Jc] = pathObj(t,x,u);

        % extract only the Jacobian w.r.t. x variables, and augment with
        % extra state assoitated with pathCost
        Jc = [Jc(1+(1:nState),:);zeros(1,length(t))];

        % reshape into the form 1 x nStateCmbn x nTime
        Jc = reshape(Jc,1,nStateCmbn,length(t));
    end

    dz = [dx;dc];

    J = cat(1,Jx,Jc);
end

end


function [F,J] = backwardEuler(z,xk,uk,uk1,dynamics,pathObj,t,dt)
% implicit integration: x[k+1] = x[k] + dt f(x[k+1],u[k+1])

z = reshape(z,size(xk,1)+1,size(xk,2));

if nargout == 1

    [fk1] = combinedDynGrad(t+dt,xk+z(1:end-1,:),uk1,dynamics,pathObj);

    F = -z + dt * (fk1);
    
    % vectorize
    F = F(:);

else

    [fk1,dfk1] = combinedJacobian(t+dt,xk+z(1:end-1,:),uk1,dynamics,pathObj);

    F = -z + dt * (fk1);
    
    % vectorize
    F = F(:);
    
    % vectorized jacobian
    dfk1 = num2cell(dfk1,[1 2]);
    J = -eye(numel(z)) + dt * blkdiag(dfk1{:});

end

end

function [F,J] = trapezoidal(z,xk,uk,uk1,dynamics,pathObj,t,dt)
% implicit integration: x[k+1] = x[k] + dt ( f(x[k+1],u[k+1]) + f([x[k],u[k]) )

z = reshape(z,size(xk,1)+1,size(xk,2));

if nargout == 1
  
    fk = combinedDynGrad(t,xk,uk,dynamics,pathObj);
    fk1 = combinedDynGrad(t+dt,xk+z(1:end-1,:),uk1,dynamics,pathObj);

    F = -z + dt/2 * (fk+fk1);
  
    % vectorize
    F = F(:);
  
else
  
    [fk] = combinedJacobian(t,xk,uk,dynamics,pathObj);
    [fk1,dfk1] = combinedJacobian(t+dt,xk+z(1:end-1,:),uk1,dynamics,pathObj);

    F = -z + dt/2 * (fk+fk1);

    % vectorize
    F = F(:);
    
    % vectorized jacobian
    dfk1 = num2cell(dfk1,[1 2]);
    J = -eye(numel(z)) + dt/2 * blkdiag(dfk1{:});
  
end
end

function [F,J] = hermiteSimpson(z,xk,uk,uMid,uk1,dynFun,pathObj,t,dt)
% hermite simpson: 
% x[k+1] = x[k] + dt/6 ( f[k] + 4fbar + f[k+1])
% fbar = f[ybar,(u[k+1]+u[k])/2)
% ybar = 0.5(y[k]+y[k+1]) + dt/8*(f[k]-f[k+1])

if nargout == 1

    fk = combinedDynGrad(t,xk,uk,dynFun,pathObj);
    fk1 = combinedDynGrad(t+dt,xk+z(1:end-1),uk1,dynFun,pathObj);

    % Something strange happening here in the combined dynamics...
    ybar = 0.5 * (z(1:end-1)+xk + xk) + dt/8 * (fk(1:end-1) - fk1(1:end-1));

    fbar = combinedDynGrad(t+0.5*dt,ybar,uMid,dynFun,pathObj);


    F = -z + dt/6 * (fk + 4*fbar + fk1);

else

    % There may be something fishy going on here. Why is the full z not
    % included in the integration?

    fk = combinedJacobian(t,xk,uk,dynFun,pathObj);
    [fk1,dfk1] = combinedJacobian(t+dt,[xk]+z(1:end-1),uk1,dynFun,pathObj);

    ybar = 0.5 * (z(1:end-1)+xk+xk) + dt/8 * (fk(1:end-1) - fk1(1:end-1));

    dybar = 0.5*eye(2,3) + dt/8 * (-dfk1(1:end-1,:));

    [fbar,dfbar] = combinedJacobian(t+0.5*dt,ybar,uMid,dynFun,pathObj);

    F = -z + dt/6 * (fk + 4*fbar + fk1);

    J = -eye(3) + dt/6 * (4*dfbar(:,1:end-1)*dybar + dfk1);

end

end

function [F,J] = hermiteSimpson_v(z,xk,uk,uMid,uk1,dynFun,pathObj,t,dt)
% hermite simpson: 
% x[k+1] = x[k] + dt/6 ( f[k] + 4fbar + f[k+1])
% fbar = f[ybar,(u[k+1]+u[k])/2)
% ybar = 0.5(y[k]+y[k+1]) + dt/8*(f[k]-f[k+1])

z = reshape(z,size(xk,1)+1,size(xk,2));

if nargout == 1

    fk = combinedDynGrad(t,xk,uk,dynFun,pathObj);
    fk1 = combinedDynGrad(t+dt,xk+z(1:end-1,:),uk1,dynFun,pathObj);

    % Something strange happening here in the combined dynamics...
    ybar = 0.5 * (z(1:end-1,:)+xk + xk) + dt/8 * (fk(1:end-1,:) - fk1(1:end-1,:));

    fbar = combinedDynGrad(t+0.5*dt,ybar,uMid,dynFun,pathObj);

    F = -z + dt/6 * (fk + 4*fbar + fk1);
    
    % vectorize
    F = F(:);

else

    % There may be something fishy going on here. Why is the full z not
    % included in the integration?

    fk = combinedJacobian(t,xk,uk,dynFun,pathObj);
    [fk1,dfk1] = combinedJacobian(t+dt,xk+z(1:end-1,:),uk1,dynFun,pathObj);

    ybar = 0.5 * (z(1:end-1,:)+xk + xk) + dt/8 * (fk(1:end-1,:) - fk1(1:end-1,:));

    dybar = 0.5*repmat(eye(2,3),1,1,length(t)) + dt/8 * (-dfk1(1:end-1,:,:));
    dybar = num2cell(dybar,[1 2]);
    dybar = blkdiag(dybar{:});

    [fbar,dfbar] = combinedJacobian(t+0.5*dt,ybar,uMid,dynFun,pathObj);

    F = -z + dt/6 * (fk + 4*fbar + fk1);
    F = F(:);

    dfbar = num2cell(dfbar(:,1:end-1,:),[1 2]);
    dfbar = blkdiag(dfbar{:});
    
    dfk1 = num2cell(dfk1,[1 2]);
    dfk1 = blkdiag(dfk1{:});
    
    J = -eye(numel(z)) + dt/6 * (4*dfbar*dybar + dfk1);

end

end
