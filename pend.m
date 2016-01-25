function result = pend(duration, targetangle, N, oldresult);

	% Finds optimal motion of a torque-driven pendulum.  Task is to move from one
	% static posture to another in a given time.
	% Method: direct collocation with 2-point discretization formula for acceleration.
	
	% Authors: Ton van den Bogert <a.vandenbogert@csuohio.edu>
    %          Milad Zarei <m.zarei@csuohio.edu>
    
    % Basically this code developed by Dr. van den Bogert, and revised by
    % Milad to caluculate velocity

	% Inputs:
	%	duration		duration of the movement (s)
	%	targetangle		target angle (radian)
	%	N				number of collocation nodes to use
	% 	oldresult		(optional) initial guess
	
	% Notes:
	% 1. This code may be useful as a template for solving other optimal control problems, such
	%    as cart-pole upswing.
	% 2. IPOPT will be used when it is installed, otherwise Newton's method.  IPOPT is recommended
	%    because it is more robust.  Newton's method can still solve most problems, but you may
	%    need to solve a sequence of problems of increasing difficulty to ensure convergence.
	% 3. The solution may be a local optimum, especially for tasks that involve multiple
	%    revolutions.  Try different initial guesses (maybe random) to check for this.
	
	% The following examples all converge with IPOPT or Newton:
	%	r = pend(1.0, pi , 100);		% swing up in 1 second, 100 collocation nodes
	%   r = pend(3.0, pi , 100);		% now do it in 3 seconds, note the countermovement
	%   r = pend(10.0, pi , 100);		% now do it in 10 seconds, multiple countermovements are seen
	%   r = pend(5.0, 2*pi, 300);		% do a full revolutions in 5 seconds, 300 collocation nodes
	%   r2 = pend(..,..,..,r);			% use previous result r as initial guess
	tic
    global model
	% settings
	MaxIterations = 1000;
	if exist('ipopt', 'file') == 3
		method = 'ipopt';
	else
		disp('IPOPT is not installed.');
		disp('Newton method will be used, may be more sensitive to initial guess.');
		disp('Hit ENTER to continue...');
		pause
		method = 'newton';
	end

	% initializations
	close all
    model.N  = N;
    model.targetangle = targetangle;
	model.h = duration/(N-1);           % time interval between nodes
	model.times = model.h*(0:N-1)';		% list of time points
	model.Ncon = 2*(N-1) + 4;           % N-1 dynamics constraints and 4 task constraints

	% model parameters
	model.L = 1;			%	length of pendulum (m)
	model.m = 1;			%	mass of pendulum (kg)
	model.I = 1;			% 	moment of inertia relative to pivot (kg m^2)
	model.g = 9.81;         %	gravity (m s^-2)
	% state variable are x & v. x: angle relative to hanging down
    % v: velocity of the pendulum (rad/sec)
	% control variable is u: torque applied at joint

	% if oldresult was provided, use it as initial guess, otherwise use zero initial guess (pendulum hangs down, no torque)
	if (nargin == 4)
        load('result')
        oldresult = result;
		oldN = numel(oldresult.t);
		oldreltime = (0:oldN-1)'/(oldN-1);                      % sample times of old result, from 0 to 1
		newreltime = (0:model.N-1)'/(model.N-1);				% sample times for new optimization, from 0 to 1
		x = interp1(oldreltime, oldresult.x, newreltime);
        v = interp1(oldreltime, oldresult.v, newreltime);
		u = interp1(oldreltime, oldresult.u, newreltime);
	else
		x = zeros(model.N,1);
        v = zeros(model.N,1);
		u = zeros(model.N,1);
	end

	% encode initial guess of unknowns into a long column vector X
	X0 = [x ; v; u];
	model.ix = (1:model.N);                   % index to elements in X where angles x are stored
	model.iv = model.N + (1:model.N);         % index to elements in X where velocities are stored
    model.iu = 2*model.N + (1:model.N);       % index to elements in X where controls u are stored
	model.NX = size(X0,1);                    % number of unknowns
	show(X0, confun(X0));                     % show the initial guess

    if N == 10
        checkderiv; % checking the derivations that function calculates with the numerical derivations
    end
	
	if strcmp(method,'ipopt')
		% solve the NLP with IPOPT
		funcs.objective         = @objfun;
		funcs.gradient          = @objgrad;
		funcs.constraints       = @confun;
		funcs.jacobian          = @conjac;
		funcs.jacobianstructure = @conjacstructure;
		options.cl              = zeros(model.Ncon,1);
		options.cu              = zeros(model.Ncon,1);	
		options.ipopt.max_iter  = MaxIterations;
		options.ipopt.hessian_approximation = 'limited-memory';
		
        
        [X, info] = ipopt(X0,funcs,options);
         comp.t=toc;
		% when user specified 1 node, we do a derivative check
     
        elseif strcmp(method, 'newton')
		solve the NLP using Newton iteration on the KKT conditions
		X = X0;
		ctol = 1e-4;		% constraint tolerance
		ftol = 1e-4;		% cost function tolerance
		xtol = 1e-4;		% solution tolerance
		F = 1e10;
		for iter=1:MaxIterations
			Fprev = F;

			% evaluate objective function F and constraint violations c
			F = objfun(X);
			G = objgrad(X);
			H = objhess(X);
			c = confun(X);
			J = conjac(X);

			% form the linearized KKT system K*x = b				
			K = [H J'; J sparse(Ncon,Ncon)];		
			b = [-G; -c];
			
			% solve the linear system K*dZ=b 
			% Z is a vector containing the unknowns X and the Lagrange multipliers.  dZ is the change in this iteration
			dZ = K\b;	
			dX = dZ(1:NX);					% the first NX are the elements of X
			
			% do a half Newton step (converges slower than full Newton step, but more likely to converge)
			% for more robust convergence, we should do a line search here to make sure we always have progress 
			X = X + dX/2;
			rmsC = sqrt(mean(c.^2));
			rmsdX = sqrt(mean(dX.^2));
			fprintf('Iter: %3d  F=%10.5e  rms(c)=%10.5e   rms(dX)=%10.5e\n', iter,F,rmsC,rmsdX);
					
			if (max(abs(c)) < ctol) && (abs(F-Fprev)<ftol) && (mean(abs(dX))<xtol)
				break;
			end 
		end
		if iter >= MaxIterations
			disp('Maximum number of iterations exceeded.');
		else
			disp('Optimal solution found');
		end
	else
		error('method not recognized');
	end
	
	% plot results
	show(X, confun(X));
	
	% make movie of the solution
	disp('Hit ENTER to generate animation...');
	pause
    avi = VideoWriter('pendDC.avi');
    avi.FrameRate = 15;
    open(avi);
	figure(2);
	clf;
	set(gcf,'Position',[5 100 650 650]);
	set(gcf, 'color', 'white');
	s = 1.5*model.L;
	for i=1:model.N
		plot([-s s],[0 0],'k','LineWidth',2);
		hold on
		plot([0 model.L*cos(X(i)-pi/2)], [0 model.L*sin(X(i)-pi/2)],'b-o','LineWidth',2);
		axis('equal');
		axis('square');
		axis([-s s -s s]);
		title(['t = ' num2str(model.times(i),'%8.3f')]);
		if (i==1)
			F = getframe(gca);
			frame = [1 1 size(F.cdata,2) size(F.cdata,1)];
		else
			F = getframe(gca,frame);
		end
		writeVideo(avi,F);
		drawnow;
		hold off;
	end
	close(avi);
    
    if info.status == 0;
        result.status = 'Solved';
        
    else info.status ~= 0;
        result.status = 'Failed';
    end

	% store results
	result.t = model.times;
	result.x = X(1:model.N);
    result.v = X(model.N+(1:model.N));
	result.u = X(2*model.N+(1:model.N));
	result.cost = model.h*(sum(result.u.^2));
    result.compt = comp.t;
    
    % save result
    filename = 'result.mat';
    save(filename,'result')
    
end

	%=========================================================
	function F = objfun(X)
	% objective function: integral of squared controls
    global model
    iu = model.iu;
	F  = sum(X(iu).^2);
	end

	%=========================================================
	function G = objgrad(X)
    % gradient of the objective function coded in objfun
    global model
    NX = model.NX;
    iu = model.iu;
    h  = model.h;
    G  = zeros(NX,1);
    G(iu) = 2 * h * X(iu);
	end

	%=========================================================
	function H = objhess(X)
    % hessian of objective function coded in objfun
    global model 
    NX = model.NX;
    N  = model.N;
    iu = model.iu;
    h  = model.h;
	H  = spalloc(NX,NX,N);
	H(iu,iu) = 2 * h * speye(N,N);
	end

	%=========================================================
	function c = confun(X)
    global model
    Ncon = model.Ncon;
    N    = model.N;
    h    = model.h;
    m    = model.m;
    L    = model.L;
    g    = model.g;
    I    = model.I;
    targetangle = model.targetangle;
		
    % constraint function (dynamics constraints and task constraints)
		
		% size of constraint vector
		c = zeros(Ncon,1);        % there are 2 state variables, so we have 2 equations. Then, the constrsints should be doubled
        
		% dynamics constraints
		% Note: torques at node 1 and node N do not affect movement and will therefore
		% always be zero in a minimal-effort solution.
		for i=1:N-1
			x1       = X(i); 
			x2       = X(i+1);
            v1       = X(N+i);
            v2       = X(N+1+i);                                         % two-point formula for angular acceleration
            u2       = X(2*N+1+i);
            c(i)     = (x2-x1)/h - v2;                                   % d(theta)/dt = thetadot
			c(N-1+i) = (v2-v1)/h - ( -m * g * L*sin(x2) + u2) / I;	     % d(thetadot)/dt = -mgLsin(theta)+u
            % the right hand side is always evaluated at second point time (x2, v2, u2) [Backward Euler, BE]
		end
		
		% task constraints
        
		% initial position must be zero:
		c(2*(N-1)+1) = X(1);					
		% initial velocity must be zero:
        c(2*(N-1)+2) = X(N+1);
		% final position must be at target angle:
		c(2*(N-1)+3) = X(N) - targetangle;	
		% final velocity must be zero:
        c(2*(N-1)+4) = X(2*N);			
	end
	%=========================================================
	function J = conjac(X)
	global model
    Ncon = model.Ncon;
    NX   = model.NX;
    N    = model.N;
    h    = model.h;
    m    = model.m;
    L    = model.L;
    g    = model.g;
    I    = model.I;
     
		% size of Jacobian
		J = spalloc(Ncon,NX,7*(N-1) + 4);

		% dynamics constraints
		for i=1:N-1
			% Jacobian matrix: derivatives of c(i) and c(N+i) with respect to the elements of X(x1, x2, v1, v2, u2)
            % x1     = X(i); 
            % x2     = X(i+1);
            % v1     = X(N+i);
            % v2     = X(N+1+i);
            % u2     = X(2*N+1+i);
            
            % c(i)   = (x2-x1)/h - v2;                               % eq(1)
            % we need to take derivative based on parameters of eq(1)    
			J(i,i) 		= -1/h;                   % derv. with respect of x1
			J(i,i+1) 	= 1/h;                    % derv. with respect of x2
			J(i,N+i+1) 	= -1;                     % derv. with respect of v2
 			
            % c(N-1+i) = (v2-v1)/h - ( -m * g * L*sin(x2) + u2) / I;   % eq(2)
            % we need to take derivative based on parameters of eq(2)
            x2               = X(i+1);
            J(N-1+i,N+i+1)   = 1/h;               % derv. with respect of v2
            J(N-1+i,N+i)     = -1/h;              % derv. with respect of v1
            J(N-1+i,i+1)     = m*g*L*cos(x2)/I;   % derv. with respect of x2
            J(N-1+i,2*N+i+1) = -1/I;              % derv. with respect of u2
		end
		
		% task constraints

		% initial position must be zero:
		J(2*(N-1)+1, 1)   = 1;
		% initial velocity must be zero:
		J(2*(N-1)+2, N+1) = 1;
		% final position must be at target angle:
		J(2*(N-1)+3, N)   = 1;
		% final velocity must be zero:
		J(2*(N-1)+4, 2*N) = 1;
	end
	%=========================================================
	function J = conjacstructure(X)
    % number of non-zeros arrayes are defined in this function
    % for instance, if there exists J(i,i), so it has to be equal 1.
    % otherwise it means that it should be zero
    global model
    Ncon = model.Ncon;
    NX   = model.NX;
    N    = model.N;

     
		% size of Jacobian
		J = spalloc(Ncon,NX,7*(N-1) + 4);

		% dynamics constraints
		for i=1:N-1
			% Jacobian matrix: derivatives of c(i) and c(N+i) with respect to the
			% elements of X(x1, x2, v1, v2, u2)
			J(i,i) 		     = 1;
			J(i,i+1) 	     = 1;
			J(i,N+i+1) 	     = 1;
            J(N-1+i,N+i+1)   = 1;
            J(N-1+i,N+i)     = 1;
            J(N-1+i,i+1)     = 1;
            J(N-1+i,2*N+i+1) = 1;
		end
		
		% task constraints

		% initial position must be zero:
		J(2*(N-1)+1, 1)   = 1;
		% initial velocity must be zero:
		J(2*(N-1)+2, N+1) = 1;
		% final position must be at target angle:
		J(2*(N-1)+3, N)   = 1;
		% final velocity must be zero:
		J(2*(N-1)+4, 2*N) = 1;
			
    end
    %================================================================
    function checkderiv
	% use finite differences to check that the code in objgrad and conjac is correct
	global model
    Ncon = model.Ncon;
    NX   = model.NX;
    N    = model.N;
    
	hh       = 1e-6;
    X        = randn(NX,1);
	f        = objfun(X);
	grad     = objgrad(X);
	c        = confun(X);
	cjac     = conjac(X);
	cjac_num = spalloc(Ncon, NX,7*(N-1) + 4);
	grad_num = zeros(NX,1);
	
    for i=1:NX
        fprintf('checking derivatives for unknown %4d of %4d\n',i,3*N);
        Xisave        = X(i);
        X(i)          = X(i) + hh;
        cjac_num(:,i) = sparse(confun(X) - c)/hh;
        grad_num(i)   = (objfun(X) - f)/hh;
        X(i)          = Xisave;
    end
	
	% report maximal differences between analytical derivatives and numerical results
	fprintf('Max. error in constraint jacobian: ');
	matcompare(cjac, cjac_num);
	fprintf('Max. error in objective gradient: ');
	matcompare(grad, grad_num);
	disp('Hit ENTER to continue');
		
    end

    %============================================================
    function matcompare(a,b)
	% compares two matrices and prints element that has greatest difference
	[maxerr,irow] = max(abs(a-b));
	[maxerr,icol] = max(maxerr);
	irow = irow(icol);
	fprintf('%9.6f at %d %d (%9.6f vs. %9.6f)\n', full(maxerr), irow, icol, full(a(irow,icol)), full(b(irow,icol)));
    end

    %============================================================
	function show(X,c)
    global model
    ix    = model.ix;
    iv    = model.iv;
    iu    = model.iu;
    times = model.times;
		% plot the current solution
		x = X(ix);
        v = X(iv);
		u = X(iu);
		figure(1)
		subplot(4,1,1);plot(times,x);title('angle')
		subplot(4,1,2);plot(times,u);title('torque');grid on;
		subplot(4,1,3);plot(times,v);title('velocity');
        subplot(4,1,4);plot(c);title('constraint violations');
    end
 
