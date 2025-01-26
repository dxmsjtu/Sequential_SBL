function [ gamma, report] = SBLseq_v3p12( A , Y, gammaIn, options )
% function [ gamma , report ] = SBL_v3p1( A , Y, options )
% The idea behind SBL is to find a diagonal replica 'covariance' Gamma. Minimizing (YY^T / AGA^T + penality) should lead to the correct
% replica selection (up to a bogus scale factor/amplitude).
% Attention: If Y is single snapshot (and single frequency), it needs to be a row vector (the code makes a 2nd snapshot with repmat).
%% Inputs
% A - Multiple frequency augmented dictionary <f , n , m>
%     f: number of frequencies
%     n: number of sensors
%     m: number of replicas
%   Note: if f==1, A = < n , m >
%% Y - Multiple snapshot multiple frequency observations <f , n , L>
%     f: number of frequencies
%     n: number of sensors
%     L: number of snapshots
%   Note: if f==1, Y = < n , L >
% options - see SBLset.m 
% Outputs
% gamma <m , 1> - vector containing source power
%                 1: surfaces found by minimum error norm
% report - various report options
%--------------------------------------------------------------------------
% Version 1.0:
% Code originally written by P. Gerstoft.
%
% Version 2.23
% Edited to include multiple frequency support: 5/16/16
%
% Version 3.1
% Different convergance norm and code update
% A and Y have now one more dimensions
% Posterior unbiased mean
% Handles single snapshot
% Version 3.12
% more efficient diagonal gamma computations
%
% Kay L Gemba & Santosh Nannuru
% MPL/SIO/UCSD gemba@ucsd.edu & snannuru@ucsd.edu

%% check function
if ismatrix(A) % SBL needs frequency dimension
    B(1,:,:) = A;
    A        = B;
end
% number of frequencies
Nfreq     = size(A,1);
% single frequency single snapshot
if ismatrix(Y)
    % either 1 freq or 1 snapshot
    if  Nfreq == 1        
        if size(Y,2) == 1 % single snapshot
            Y=Y.'; %squeeze2
        else
            Y = permute(Y,[ 3 1 2 ]); % works
        end        
    end    
end
%%
options.SBL_v = '3.12';
%% slicing
Nsource = options.Nsource;
if options.tic == 1
    %tic
end
%% Initialize variables
% number of sensors
Nsensor   = size(A,2);
% number of dictionary entries
Ntheta    = size(A,3);
% number of snapshots in the data covariance
Nsnapshot = size(Y,3);
% noise power initialization
sigc      = ones(Nfreq,1) * options.noisepower.guess;
% posterior
x_post    = zeros(Nfreq, Ntheta, Nsnapshot);
% minimum (global) gamma
gmin_global   = realmax;
% L1 error
errornorm    = zeros(options.convergence.maxiter,1);

% initialize equal and uncorrelated weights
% gamma        = 1*ones(Ntheta,1);
gamma_num    = zeros(Nfreq , Ntheta);
gamma_denum  = zeros(Nfreq , Ntheta);

% Sample Covariance Matrix
SCM = zeros( Nfreq , Nsensor , Nsensor );

for i_f = 1 : Nfreq
    SCM(i_f,:,:) = squeeze2(Y(i_f,:,:)) * squeeze2(Y(i_f,:,:))' / Nsnapshot;
end
%% Initialize for Sequential processing change by YP
% Simple source motion model (Random walk model with finite lags alphalag)
alphalags = 2;  % cover [-alphalag , +alphalag] grid bin walk
                % ex. Beamforming with 0.5[deg] angular grid bin with lag 2
                % the previous estimate 0 deg -> informative to the next [-1:0.5:1] [deg]

% alpha distribution
alphalength = 2*alphalags+1; alphatmp = ones(1,alphalength);  %uniform
% alphalength = 2*alphalags+1; alphatmp = rand(1,alphalength);  %random
% alphatmp = normpdf(-alphalags:alphalags,0,1);                 %gaussian, here with variance 1
alpha = alphatmp / sum(alphatmp);

Ftmp = zeros(Ntheta,1); Ftmp(1:alphalags+1) = alpha(alphalags+1:2*alphalags+1);

F = toeplitz(Ftmp);
% F = eye(Ntheta);

% gammaIn = conv(gammaIn, power(alpha,2)); gammaIn(1:alphalags)=[]; gammaIn(end-alphalags+1:end)=[];
gammaIn = diag(F * diag(gammaIn) * F.');

gamma        = 1*ones(Ntheta,1) * max(gammaIn)/1000;
% Af = squeeze(A(:,:,iF));
rho = 1.0;
gamma = rho*rho*gammaIn + gamma;

%% Main Loop
%display(['Sequential SBL version ', options.SBL_v ,' initialized.']);

for j1 = 1 : options.convergence.maxiter    
    % for error analysis
    gammaOld = gamma;    
    %% gamma update    
    for i_f = 1 : Nfreq        
        Af = squeeze(A(i_f,:,:));        
        ApSigmaYinv  = Af' / (sigc(i_f) * eye(Nsensor) + Af * (repmat(gamma, [1 Nsensor] ) .* Af'));        
        % Sum over snapshots and normalize, abs for roundoff errors
        gamma_num(i_f,:)   = sum ( abs ( ( ApSigmaYinv * squeeze2(Y(i_f,:,:)) ).^2 ),2 ) / Nsnapshot;        
        % positive def quantity, abs for roundoff errors
        gamma_denum(i_f,:) = abs( sum  ( ApSigmaYinv.' .* Af, 1 ) );
    end    
    % Fixed point Eq. update
    gamma = gamma  .* ((sum( gamma_num   ,1 ) ./sum( gamma_denum ,1 ) ).^(1/options.fixedpoint) ).' ;    
    %% sigma and L2 error using unbiased posterior update     
    % locate same peaks for all frequencies
    [ ~ , Ilocs] = findpeaks(gamma,'SORTSTR','descend','NPEAKS',Nsource);
    Apeak      = A(:,:,Ilocs);    
    for i_f = 1 : Nfreq        
        % only active replicas
        Am     = squeeze2(Apeak(i_f,:,:));   
        % noise estimate
        sigc(i_f) = real(trace( (eye(Nsensor)-Am*pinv(Am)) * squeeze(SCM(i_f,:,:)) ) / ( Nsensor- Nsource ) );
    end    
    %% Convergance
    % checks convergance and displays status reports    
    % convergance indicator
    errornorm(j1) = norm ( gamma - gammaOld, 1 ) / norm ( gamma, 1 );
    
    % look into the past and find best error since then
    if j1 > options.convergence.min_iteration  &&  errornorm(j1) < gmin_global
        gmin_global = errornorm(j1);
        gamma_min   = gamma;
        iteration_L1   = j1;
    end
    
    % inline convergence code
    if j1 > options.convergence.min_iteration && ( errornorm(j1) < options.convergence.error  || iteration_L1 + options.convergence.delay <= j1)
        
        if options.flag == 1
            display(['Solution converged. Iteration: ',num2str(sprintf('%.4u',j1)),'. Error: ',num2str(sprintf('%1.2e' , errornorm(j1) )),'.'])
        end
        break; % goodbye
        
        % not convereged
    elseif j1 == options.convergence.maxiter
        warning(['Solution not converged. Error: ',num2str(sprintf('%1.2e' , errornorm(j1) )),'.'])        
        % status report
    elseif j1 ~= options.convergence.maxiter  && options.flag == 1 && mod(j1,options.status_report) == 0 % Iteration reporting
        display(['Iteration: ',num2str(sprintf('%.4u',j1)),'. Error: ',num2str(sprintf('%1.2e' , errornorm(j1) )),'.' ])        
    end    
end
%% Posterior distribution
% x_post - posterior unbiased mean
for i_f = 1 : Nfreq    
    Af = squeeze(A(i_f,:,:));    
    % c.f. (31) of [2]
    x_post(i_f,:,:) = repmat(gamma, [1 Nsnapshot] ) .* (Af' / (sigc(i_f) * eye(Nsensor) + Af * (repmat(gamma, [1 Nsensor] ) .* Af')) * squeeze2(Y(i_f,:,:)));    
end

%% function return
% Globla minimum
gamma = gamma_min;

% Report section

% vectors containing errors
report.results.error    = errornorm;

% Error when minimum was obtained
report.results.iteration_L1 = iteration_L1;

% General info
report.results.final_iteration.iteration = j1;
report.results.final_iteration.noisepower = sigc;

if options.tic == 1
    report.results.toc = toc;
else
    report.results.toc = 0;
end

% data
report.results.final_iteration.gamma  = gamma  ;
report.results.final_iteration.x_post = x_post ;

report.options = options;

end

function b = squeeze2(a)
%   SQUEEZE2
%   Just as squeeze but with transpose to accomodate single snapshot case.
%   This is required because matlab does not allow a singleton dimension
%   at the 'end', e.g., 3x5x1.
%   The entire code might alternatively be re-written to have the number
%   of sensors n~=1 at the end (making this fix obsolete).

if ~ismatrix(a)
    siz = size(a);
    siz(siz==1) = []; % Remove singleton dimensions.
    siz = [siz ones(1,2-length(siz))]; % Make sure siz is at least 2-D
    b = reshape(a,siz);
else
    b = a.';
end

end

%EOF