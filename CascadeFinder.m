clc; clear; close all;

%% ========================================================================
% MANUAL CASCADE PID TUNER
% SEARCH ALGORITHM: RESTARTED CMA-ES + FINAL NELDER-MEAD
%
% The Simulink model and controller equations are unchanged:
%
%   Solver type     = Fixed-step
%   Solver          = ode4 (Runge-Kutta)
%   Fixed-step size = 0.001 s
%
%   Outer PID:
%       ex = x_ref - x
%       uOuter = Kpx*ex + Kix*Integral(ex) + Kdx*Derivative(ex)
%       theta_ref = -Saturation(uOuter, +/-0.3)
%
%   Inner PID:
%       eTheta = theta_ref - theta
%       F = Saturation(Kpp*eTheta + Kip*Integral(eTheta) ...
%                      + Kdp*Derivative(eTheta), +/-20)
%
% This version removes the old Differential Evolution search, pole gates,
% artificial gain bounds, restart fractions and weighted hard walls.
%
% Gains are searched in unbounded logarithmic coordinates:
%
%       K = exp(z)
%
% Only a very wide numerical guard on z is used to prevent exp overflow.
%
% CMA-ES automatically learns:
%   - search step size
%   - gain correlations
%   - useful search directions
%   - different scales of the six PID gains
%
% A candidate is ranked by:
%   1) total normalized requirement violation
%   2) transient quality after it becomes feasible
%
% The pole calculation is retained only as a diagnostic and a smooth
% instability term. It is NOT used as a hard gate.
%
% Gain order:
%   K = [Kpx Kix Kdx Kpp Kip Kdp]
%% ========================================================================

%% =============================== CONFIG ================================

% Linear CTMS model matching the supplied P_cart and P_pend
cfg.M = 0.5;
cfg.m = 0.2;
cfg.b = 0.1;
cfg.I = 0.006;
cfg.g = 9.8;
cfg.l = 0.3;

cfg.q = (cfg.M + cfg.m)*(cfg.I + cfg.m*cfg.l^2) ...
      - (cfg.m*cfg.l)^2;

% Simulink-equivalent fixed-step ode4
cfg.dt = 0.001;

% Two-stage search
cfg.quickEnd = 4.0;
cfg.finalEnd = 10.0;

% Long validation is used only after tuning, so a delayed unstable mode
% cannot hide behind a good-looking ten-second response.
cfg.validationEnd = 30.0;

% Robust continuous-time pole margin.
% The previous code accepted poles such as -3e-9, which are effectively
% marginal and can flip unstable in Simulink because of numerical and
% realization differences.
cfg.robustPoleMargin = 0.25;

% Constant position reference
cfg.xRef = 0.1;

% Initial conditions of Transfer Fcn and Integrator blocks
cfg.x0 = 0;
cfg.xd0 = 0;
cfg.theta0 = 0;
cfg.thetad0 = 0;
cfg.Ix0 = 0;
cfg.Itheta0 = 0;

% Saturations
cfg.thetaRefSat = 0.80;
cfg.forceSat = 20;

% Required performance
cfg.TsReq = 3.0;
cfg.xTolerance = 0.02*abs(cfg.xRef);
cfg.thetaTolerance = 0.01;
cfg.xEssReq = 0.002;
cfg.thetaLimit = 20*pi/180;

cfg.maxOvershootPercent = 5;
cfg.maxOvershoot = ...
    cfg.maxOvershootPercent/100*abs(cfg.xRef);

% CMA-ES stage 1: fast exploration
cfg.quick.lambda = 18;
cfg.quick.maxGenerations = 20;
cfg.quick.sigma0 = 1.20;
cfg.quick.noImproveStop = 14;

% CMA-ES stage 2: full ten-second optimization
cfg.full.lambda = 24;
cfg.full.maxGenerations = 30;
cfg.full.sigma0 = 0.65;
cfg.full.noImproveStop = 16;
cfg.full.restarts = 3;

% Final built-in Nelder-Mead refinement
cfg.useFinalFminsearch = false;
cfg.fminMaxIter = 450;
cfg.fminMaxFunEvals = 2500;

% Numerical protection only, not an engineering gain bound:
% exp(-25) ~= 1.39e-11 and exp(25) ~= 7.20e10
cfg.maxAbsLogGain = 25;

% Stop a simulation early after obvious divergence
cfg.maxAbsX = 5;
cfg.maxAbsXd = 100;
cfg.maxAbsTheta = 2;
cfg.maxAbsThetaDot = 200;
cfg.maxAbsIntegral = 200;

rng('shuffle');

% Starting points only. They are not accepted without simulation.
seedK = [
    0.748861745  5.06526222978174e-7  0.208762347 ...
    11.4473746    2.5238694e-5         1.97890723;

    1.28281046    1.45779651           0.0524047234 ...
    74.0151638    215.804351           5.93610049;

    0.00280781562797  0.000161431408239  0.0183171639775 ...
    19.5358168162     29.1726351446      1.99833742926;

    0.10  0.01   0.05   20   2    2;
    0.50  0.10   0.15   30  10    3;
    1.00  0.50   0.10   60  50    5
];

fprintf('\n============================================================\n');
fprintf('MANUAL CASCADE PID - CMA-ES ODE4 TUNER\n');
fprintf('Search algorithm  = restarted CMA-ES\n');
fprintf('Final refinement  = Nelder-Mead fminsearch\n');
fprintf('Gain bounds       = NONE, logarithmic search\n');
fprintf('Pole hard gate    = NONE\n');
fprintf('Controller        = manual P + I + D\n');
fprintf('Derivative filter = NONE\n');
fprintf('N                 = NONE\n');
fprintf('Solver            = fixed-step ode4 / RK4\n');
fprintf('Fixed step        = %.6g s\n',cfg.dt);
fprintf('Quick horizon     = %.6g s\n',cfg.quickEnd);
fprintf('Final horizon     = %.6g s\n',cfg.finalEnd);
fprintf('Validation horizon= %.6g s\n',cfg.validationEnd);
fprintf('Pole margin       = max Re(p) <= -%.6g\n', ...
    cfg.robustPoleMargin);
fprintf('theta_ref limit   = +/- %.6g rad\n',cfg.thetaRefSat);
fprintf('Force limit       = +/- %.6g N\n',cfg.forceSat);
fprintf('Overshoot limit   = %.6g %%\n',cfg.maxOvershootPercent);
fprintf('============================================================\n\n');
drawnow;

%% =========================== CHOOSE START ===============================

fprintf('Evaluating starting points on the quick horizon...\n');
drawnow;

bestQuick = emptyResult(cfg,cfg.quickEnd);

for i = 1:size(seedK,1)

    K = seedK(i,:);

    [fitness,met,sim] = ...
        evaluateCandidate(K,cfg,cfg.quickEnd);

    fprintf(['Seed %d/%d | feasible=%d | violation=%.6g | ', ...
             'fitness=%.8g | Ts=%.5g | OS=%.5g%% | ', ...
             'ess=%.5g | maxRePole=%.5g\n'], ...
             i,size(seedK,1),met.pass,met.violation,fitness, ...
             met.TsOverall,met.OSxPercent,met.essX, ...
             met.maxPoleReal);
    drawnow;

    if fitness < bestQuick.fitness
        bestQuick.K = K;
        bestQuick.fitness = fitness;
        bestQuick.met = met;
        bestQuick.sim = sim;
    end
end

fprintf('\nBest initial point:\n');
printResult(bestQuick);
drawnow;

%% ======================== QUICK CMA-ES STAGE ============================

fprintf('\nStarting quick CMA-ES search...\n');
drawnow;

quickMean = log(bestQuick.K(:));

bestQuick = runCMAES( ...
    quickMean, ...
    cfg.quick.sigma0, ...
    cfg.quick.lambda, ...
    cfg.quick.maxGenerations, ...
    cfg.quick.noImproveStop, ...
    cfg.quickEnd, ...
    cfg, ...
    bestQuick, ...
    'QUICK');

%% ===================== FULL-HORIZON VALIDATION ==========================

fprintf('\nValidating the quick-stage controller over %.3f seconds...\n', ...
    cfg.finalEnd);
drawnow;

[fullFitness,fullMet,fullSim] = ...
    evaluateCandidate(bestQuick.K,cfg,cfg.finalEnd);

bestFull.K = bestQuick.K;
bestFull.fitness = fullFitness;
bestFull.met = fullMet;
bestFull.sim = fullSim;

fprintf('\nFull-horizon result before full search:\n');
printResult(bestFull);
drawnow;

%% ====================== RESTARTED FULL CMA-ES ===========================

for restart = 1:cfg.full.restarts

    if restart==1
        sigmaStart = cfg.full.sigma0;
    else
        % Broaden later restarts if the first basin was unhelpful.
        sigmaStart = min( ...
            1.6, ...
            cfg.full.sigma0*1.55^(restart-1));
    end

    lambda = round( ...
        cfg.full.lambda*1.35^(restart-1));

    fprintf(['\nStarting FULL CMA-ES restart %d/%d | ', ...
             'lambda=%d | sigma=%.5g\n'], ...
             restart,cfg.full.restarts,lambda,sigmaStart);
    drawnow;

    % Restart around the best point found so far, plus a small random
    % displacement after the first run.
    meanZ = log(bestFull.K(:));

    if restart>1
        meanZ = meanZ ...
              + 0.35*randn(size(meanZ));
    end

    meanZ = guardLogVector(meanZ,cfg.maxAbsLogGain);

    bestFull = runCMAES( ...
        meanZ, ...
        sigmaStart, ...
        lambda, ...
        cfg.full.maxGenerations, ...
        cfg.full.noImproveStop, ...
        cfg.finalEnd, ...
        cfg, ...
        bestFull, ...
        sprintf('FULL-R%d',restart));

    if bestFull.met.pass
        fprintf(['A time-domain and pole-robust controller has been ', ...
                 'found.\n']);
        drawnow;
        break;
    end
end

%% ======================= FINAL LOCAL REFINEMENT =========================

if cfg.useFinalFminsearch

    fprintf('\nStarting final Nelder-Mead refinement...\n');
    drawnow;

    z0 = log(bestFull.K(:)).';

    objective = @(z) objectiveOnly( ...
        z,cfg,cfg.finalEnd);

    options = optimset( ...
        'Display','iter', ...
        'MaxIter',cfg.fminMaxIter, ...
        'MaxFunEvals',cfg.fminMaxFunEvals, ...
        'TolX',1e-6, ...
        'TolFun',1e-6);

    [zCandidate,~] = ...
        fminsearch(objective,z0,options);

    zCandidate = guardLogVector( ...
        zCandidate(:),cfg.maxAbsLogGain);

    KCandidate = exp(zCandidate).';

    [candidateFitness,candidateMet,candidateSim] = ...
        evaluateCandidate(KCandidate,cfg,cfg.finalEnd);

    if candidateFitness < bestFull.fitness
        bestFull.K = KCandidate;
        bestFull.fitness = candidateFitness;
        bestFull.met = candidateMet;
        bestFull.sim = candidateSim;
    end
end


%% ======================== LONG-HORIZON VALIDATION ========================

fprintf('\nValidating the selected controller over %.3f seconds...\n', ...
    cfg.validationEnd);
drawnow;

[validationFitness,validationMet,validationSim] = ...
    evaluateCandidate(bestFull.K,cfg,cfg.validationEnd);

validationPass = validationMet.pass;

fprintf('\nLong-horizon validation result:\n');
fprintf('Validation pass       = %d\n',validationPass);
fprintf('Validation Ts overall = %.10g s\n',validationMet.TsOverall);
fprintf('Validation overshoot  = %.10g %%\n',validationMet.OSxPercent);
fprintf('Validation ess        = %.12g m\n',validationMet.essX);
fprintf('Validation max Re pole= %.12g 1/s\n', ...
    validationMet.maxPoleReal);
drawnow;

% If the 30-second test fails, continue CMA-ES on the long horizon instead
% of exporting a controller that only looks good briefly.
if ~validationPass

    fprintf(['\nTen-second tuning passed but long validation failed. ', ...
             'Starting one long-horizon CMA-ES repair run...\n']);
    drawnow;

    bestValidation.K = bestFull.K;
    bestValidation.fitness = validationFitness;
    bestValidation.met = validationMet;
    bestValidation.sim = validationSim;

    bestValidation = runCMAES( ...
        log(bestFull.K(:)), ...
        0.40, ...
        28, ...
        60, ...
        14, ...
        cfg.validationEnd, ...
        cfg, ...
        bestValidation, ...
        'VALIDATION');

    bestFull = bestValidation;

    [validationFitness,validationMet,validationSim] = ...
        evaluateCandidate(bestFull.K,cfg,cfg.validationEnd);

    validationPass = validationMet.pass;
end

%% ============================== FINAL OUTPUT ============================

fprintf('\n============================================================\n');

if bestFull.met.pass && validationPass
    fprintf('FINAL ROBUST CASCADE PID GAINS\n');
else
    fprintf('BEST CONTROLLER FOUND - ROBUST REQUIREMENTS NOT ALL MET\n');
end

fprintf('============================================================\n');

Kbest = bestFull.K;
finalFitness = bestFull.fitness;
finalMet = bestFull.met;
finalSim = bestFull.sim;

% Re-evaluate the final gains on the standard ten-second horizon for the
% main plots, while keeping the separate 30-second validation result.
[plotFitness,plotMet,plotSim] = ...
    evaluateCandidate(Kbest,cfg,cfg.finalEnd);

if all(isfinite(plotSim.x))
    finalSim = plotSim;
    finalFitness = plotFitness;
    finalMet = plotMet;
end

fprintf('\nOuter manual PID:\n');
fprintf('Kpx = %.15g\n',Kbest(1));
fprintf('Kix = %.15g\n',Kbest(2));
fprintf('Kdx = %.15g\n',Kbest(3));

fprintf('\nInner manual PID:\n');
fprintf('Kpp = %.15g\n',Kbest(4));
fprintf('Kip = %.15g\n',Kbest(5));
fprintf('Kdp = %.15g\n',Kbest(6));

fprintf('\nPerformance:\n');
fprintf('Ts_x                  = %.10g s\n',finalMet.TsX);
fprintf('Ts_theta              = %.10g s\n',finalMet.TsTheta);
fprintf('Ts_overall            = %.10g s\n',finalMet.TsOverall);
fprintf('Cart steady-state err = %.12g m\n',finalMet.essX);
fprintf('Cart overshoot        = %.12g m\n',finalMet.OSx);
fprintf('Cart overshoot        = %.12g %%\n',finalMet.OSxPercent);
fprintf('Maximum |theta|       = %.12g rad\n',finalMet.thetaMax);
fprintf('Maximum |theta_ref|   = %.12g rad\n',finalMet.thetaRefMax);
fprintf('Maximum |F|           = %.12g N\n',finalMet.Fmax);
fprintf('Force saturation      = %.8g %%\n', ...
    100*finalMet.forceSatRatio);
fprintf('ThetaRef saturation   = %.8g %%\n', ...
    100*finalMet.thetaRefSatRatio);
fprintf('Maximum real pole     = %.12g 1/s\n',finalMet.maxPoleReal);
fprintf('Required pole margin  = -%.12g 1/s\n',cfg.robustPoleMargin);
fprintf('Pole margin satisfied = %d\n',finalMet.poleMarginPass);
fprintf('30 s validation pass  = %d\n',validationPass);
fprintf('Normalized violation  = %.12g\n',finalMet.violation);
fprintf('Final fitness         = %.12g\n',finalFitness);
fprintf('All requirements met  = %d\n', ...
    finalMet.pass && validationPass);

fprintf('\nClosed-loop poles:\n');
disp(finalMet.poles.');

fprintf('============================================================\n');

%% =============================== EXPORT =================================

Kpx = Kbest(1);
Kix = Kbest(2);
Kdx = Kbest(3);

Kpp = Kbest(4);
Kip = Kbest(5);
Kdp = Kbest(6);

outerDirection = -1;

assignin('base','Kpx',Kpx);
assignin('base','Kix',Kix);
assignin('base','Kdx',Kdx);
assignin('base','Kpp',Kpp);
assignin('base','Kip',Kip);
assignin('base','Kdp',Kdp);
assignin('base','outerDirection',outerDirection);

save('Manual_Cascade_PID_CMAES_Robust_ODE4_Final.mat', ...
    'Kpx','Kix','Kdx','Kpp','Kip','Kdp', ...
    'outerDirection','cfg','finalMet','finalFitness', ...
    'validationMet','validationFitness','validationPass');

%% ================================ PLOTS =================================

figure('Name','Manual Cascade PID - CMA-ES ode4','Color','w');
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nexttile;
plot(finalSim.t,finalSim.x,'LineWidth',1.4); hold on;
plot(finalSim.t,cfg.xRef*ones(size(finalSim.t)), ...
    '--','LineWidth',1.2);
yline(cfg.xRef+cfg.maxOvershoot,'--');
grid on;
xlabel('Time (s)');
ylabel('x (m)');
title('Cart Position');
legend('x','x_{ref}','overshoot limit','Location','best');

nexttile;
plot(finalSim.t,finalSim.theta,'LineWidth',1.4); hold on;
yline(cfg.thetaLimit,'--');
yline(-cfg.thetaLimit,'--');
yline(cfg.thetaTolerance,':');
yline(-cfg.thetaTolerance,':');
grid on;
xlabel('Time (s)');
ylabel('\theta (rad)');
title('Pendulum Angle');

nexttile;
plot(finalSim.t,finalSim.thetaRef,'LineWidth',1.4); hold on;
yline(cfg.thetaRefSat,'--');
yline(-cfg.thetaRefSat,'--');
grid on;
xlabel('Time (s)');
ylabel('\theta_{ref} (rad)');
title('Outer PID Output');

nexttile;
plot(finalSim.t,finalSim.F,'LineWidth',1.4); hold on;
yline(cfg.forceSat,'--');
yline(-cfg.forceSat,'--');
grid on;
xlabel('Time (s)');
ylabel('F (N)');
title('Control Force');

sgtitle('Manual Cascade PID, CMA-ES, Fixed-Step ode4');

%% ========================================================================
%                               FUNCTIONS
% ========================================================================

function best = runCMAES( ...
    meanZ,sigma,lambda,maxGenerations,noImproveStop, ...
    horizon,cfg,best,label)

    n = numel(meanZ);

    meanZ = guardLogVector( ...
        meanZ(:),cfg.maxAbsLogGain);

    mu = floor(lambda/2);

    weights = log(mu+0.5)-log(1:mu);
    weights = weights(:)/sum(weights);

    muEff = 1/sum(weights.^2);

    cc = (4+muEff/n)/(n+4+2*muEff/n);
    cs = (muEff+2)/(n+muEff+5);

    c1 = 2/((n+1.3)^2+muEff);

    cmu = min( ...
        1-c1, ...
        2*(muEff-2+1/muEff)/((n+2)^2+muEff));

    damping = ...
        1 ...
      + 2*max(0,sqrt((muEff-1)/(n+1))-1) ...
      + cs;

    pc = zeros(n,1);
    ps = zeros(n,1);

    C = eye(n);
    B = eye(n);
    D = ones(n,1);
    invSqrtC = eye(n);

    chiN = sqrt(n)*(1-1/(4*n)+1/(21*n^2));

    noImprove = 0;

    for generation = 1:maxGenerations

        arz = randn(n,lambda);
        ary = B*(D.*arz);

        arx = meanZ+sigma*ary;

        for k = 1:lambda
            arx(:,k) = guardLogVector( ...
                arx(:,k),cfg.maxAbsLogGain);
        end

        fitness = inf(lambda,1);
        metrics = cell(lambda,1);
        simulations = cell(lambda,1);

        for k = 1:lambda

            K = exp(arx(:,k)).';

            [fitness(k),metrics{k},simulations{k}] = ...
                evaluateCandidate(K,cfg,horizon);
        end

        [fitness,order] = sort(fitness,'ascend');

        arx = arx(:,order);
        ary = ary(:,order);
        arz = arz(:,order);

        metrics = metrics(order);
        simulations = simulations(order);

        previousMean = meanZ;

        meanZ = arx(:,1:mu)*weights;
        meanZ = guardLogVector( ...
            meanZ,cfg.maxAbsLogGain);

        yMean = (meanZ-previousMean)/sigma;

        ps = (1-cs)*ps ...
           + sqrt(cs*(2-cs)*muEff) ...
           * invSqrtC*yMean;

        hsigDenominator = sqrt( ...
            1-(1-cs)^(2*generation));

        hsig = ...
            norm(ps)/max(hsigDenominator,eps)/chiN ...
            < (1.4+2/(n+1));

        pc = (1-cc)*pc ...
           + hsig*sqrt(cc*(2-cc)*muEff) ...
           * yMean;

        selectedSteps = ...
            (arx(:,1:mu)-previousMean)/sigma;

        C = ...
            (1-c1-cmu)*C ...
          + c1*( ...
                pc*pc.' ...
              + (1-hsig)*cc*(2-cc)*C) ...
          + cmu*selectedSteps ...
            * diag(weights) ...
            * selectedSteps.';

        C = (C+C.')/2;

        sigma = sigma*exp( ...
            (cs/damping)*(norm(ps)/chiN-1));

        % Numerical covariance repair
        [B,eigenValues] = eig(C);
        eigenValues = real(diag(eigenValues));
        eigenValues = max(eigenValues,1e-14);

        [eigenValues,sortIndex] = ...
            sort(eigenValues,'ascend');

        B = real(B(:,sortIndex));
        D = sqrt(eigenValues);

        C = B*diag(D.^2)*B.';
        C = (C+C.')/2;

        invSqrtC = ...
            B*diag(1./D)*B.';

        generationBestFitness = fitness(1);
        generationBestMet = metrics{1};
        generationBestSim = simulations{1};
        generationBestK = exp(arx(:,1)).';

        if generationBestFitness < best.fitness

            best.K = generationBestK;
            best.fitness = generationBestFitness;
            best.met = generationBestMet;
            best.sim = generationBestSim;

            noImprove = 0;
        else
            noImprove = noImprove+1;
        end

        fprintf(['%s gen %3d/%3d | feasible=%d | ', ...
                 'violation=%.6g | fitness=%.8g | ', ...
                 'Ts=%.5g | OS=%.5g%% | ess=%.5g | ', ...
                 'maxRePole=%.5g | poleOK=%d | ', ...
                 'sigma=%.5g | noImprove=%d\n'], ...
                 label,generation,maxGenerations, ...
                 best.met.pass,best.met.violation, ...
                 best.fitness,best.met.TsOverall, ...
                 best.met.OSxPercent,best.met.essX, ...
                 best.met.maxPoleReal,best.met.poleMarginPass, ...
                 sigma,noImprove);
        drawnow;

        if best.met.pass && noImprove>=noImproveStop
            break;
        end

        if noImprove>=2*noImproveStop && sigma<0.03
            break;
        end
    end
end

function value = objectiveOnly(z,cfg,horizon)

    z = guardLogVector( ...
        z(:),cfg.maxAbsLogGain);

    K = exp(z).';

    value = evaluateCandidate(K,cfg,horizon);
end

function [fitness,met,sim] = ...
    evaluateCandidate(K,cfg,horizon)

    if any(~isfinite(K)) || any(K<=0)

        fitness = 1e15;
        met = invalidMetrics();
        sim = invalidSimulation(cfg,horizon);
        return;
    end

    % Pole information is diagnostic and a smooth instability indicator.
    [~,poles,maxPoleReal,charPoly] = ...
        closedLoopPoleCheck(K,cfg);

    sim = simulateExactODE4(K,cfg,horizon);
    [met,ok] = calculateMetrics(sim,cfg,horizon);

    met.poles = poles;
    met.maxPoleReal = maxPoleReal;
    met.characteristicPolynomial = charPoly;

    if isfinite(maxPoleReal)
        met.poleMarginViolation = max( ...
            0, ...
            (maxPoleReal+cfg.robustPoleMargin) ...
            / cfg.robustPoleMargin);
    else
        met.poleMarginViolation = inf;
    end

    met.poleMarginPass = ...
        isfinite(maxPoleReal) && ...
        maxPoleReal<=-cfg.robustPoleMargin;

    if ~ok

        poleViolation = max(0,maxPoleReal);

        met.poleMarginPass = false;
        met.poleMarginViolation = inf;

        met.violation = ...
            100+10*poleViolation^2;

        met.quality = 1e6;

        fitness = ...
            1e8+1e5*met.violation;

        return;
    end

    late = sim.t>=min(cfg.TsReq,horizon);

    maxLateX = max( ...
        abs(sim.x(late)-cfg.xRef));

    maxLateTheta = max( ...
        abs(sim.theta(late)));

    vSettleX = max( ...
        0,maxLateX/cfg.xTolerance-1);

    vSettleTheta = max( ...
        0,maxLateTheta/cfg.thetaTolerance-1);

    vOvershoot = max( ...
        0,met.OSx/cfg.maxOvershoot-1);

    vEss = max( ...
        0,met.essX/cfg.xEssReq-1);

    vTheta = max( ...
        0,met.thetaMax/cfg.thetaLimit-1);

    vForce = max( ...
        0,met.Fmax/cfg.forceSat-1);

    % Robust pole violation. This is not an arbitrary speed target:
    % it prevents effectively marginal solutions such as maxRe=-3e-9,
    % which were the exact cause of the delayed Simulink divergence.
    vPoleMargin = met.poleMarginViolation;

    met.violation = ...
        60*vSettleX^2 ...
      + 20*vSettleTheta^2 ...
      + 35*vOvershoot^2 ...
      + 25*vEss^2 ...
      + 8*vTheta^2 ...
      + 5*vForce^2 ...
      + 80*vPoleMargin^2;

    % A controller is feasible only if both the time response and the
    % robust pole margin pass.
    met.pass = met.pass && met.poleMarginPass;

    forceRMS = sqrt(mean(sim.F.^2));

    if isfinite(met.TsOverall)
        TsQuality = met.TsOverall;
    else
        TsQuality = ...
            horizon ...
          + maxLateX/max(cfg.xTolerance,eps) ...
          + maxLateTheta/max(cfg.thetaTolerance,eps);
    end

    met.quality = ...
        TsQuality ...
      + 0.12*met.OSxPercent ...
      + 2.0*met.essX/max(cfg.xEssReq,eps) ...
      + 0.015*forceRMS ...
      + 0.20*met.forceSatRatio ...
      + 0.10*met.thetaRefSatRatio;

    % Feasibility-first scalar ranking for CMA-ES.
    if met.pass
        fitness = met.quality;
    else
        fitness = ...
            1e5 ...
          + 1e4*met.violation ...
          + met.quality;
    end
end

function [isStable,poles,maxReal,charPoly] = ...
    closedLoopPoleCheck(K,cfg)

    Kpx = K(1);
    Kix = K(2);
    Kdx = K(3);

    Kpp = K(4);
    Kip = K(5);
    Kdp = K(6);

    q = cfg.q;

    a3 = cfg.b*(cfg.I+cfg.m*cfg.l^2)/q;
    a2 = -(cfg.M+cfg.m)*cfg.m*cfg.g*cfg.l/q;
    a1 = -cfg.b*cfg.m*cfg.g*cfg.l/q;

    Dplant = [1 a3 a2 a1 0];

    Ncart = [ ...
        (cfg.I+cfg.m*cfg.l^2)/q, ...
        0, ...
        -cfg.m*cfg.g*cfg.l/q];

    Npend = [cfg.m*cfg.l/q 0 0];

    Nx = [Kdx Kpx Kix];
    Ntheta = [Kdp Kpp Kip];

    term1 = conv([1 0 0],Dplant);
    term2 = conv(conv(Ntheta,Nx),Ncart);
    term3 = conv(conv([1 0],Ntheta),Npend);

    charPoly = addAlignedPolynomials( ...
        term1,-term2,term3);

    scale = max(abs(charPoly));

    if isempty(scale) || ...
       ~isfinite(scale) || scale==0

        poles = NaN(6,1);
        maxReal = inf;
        isStable = false;
        return;
    end

    firstValid = find( ...
        abs(charPoly)>max(1e-12*scale,eps), ...
        1,'first');

    if isempty(firstValid)

        poles = NaN(6,1);
        maxReal = inf;
        isStable = false;
        return;
    end

    charPoly = charPoly(firstValid:end);

    if numel(charPoly)~=7 || ...
       ~all(isfinite(charPoly)) || ...
       abs(charPoly(1))<1e-12

        poles = NaN(6,1);
        maxReal = inf;
        isStable = false;
        return;
    end

    poles = roots(charPoly);

    if isempty(poles) || ...
       ~all(isfinite(poles))

        maxReal = inf;
        isStable = false;
        return;
    end

    maxReal = max(real(poles));
    isStable = maxReal<0;
end

function p = addAlignedPolynomials(varargin)

    maximumLength = 0;

    for i = 1:nargin
        maximumLength = max( ...
            maximumLength,numel(varargin{i}));
    end

    p = zeros(1,maximumLength);

    for i = 1:nargin

        term = varargin{i};
        first = maximumLength-numel(term)+1;

        p(first:end) = ...
            p(first:end)+term;
    end
end

function sim = simulateExactODE4(K,cfg,horizon)

    dt = cfg.dt;
    t = (0:dt:horizon).';
    n = numel(t);

    X = zeros(n,6);

    X(1,:) = [ ...
        cfg.x0, ...
        cfg.xd0, ...
        cfg.theta0, ...
        cfg.thetad0, ...
        cfg.Ix0, ...
        cfg.Itheta0];

    thetaRef = zeros(n,1);
    force = zeros(n,1);

    valid = true;

    previousMajorEx = 0;
    previousMajorETheta = 0;

    for k = 1:n-1

        tk = t(k);
        Xk = X(k,:).';

        initialMajor = (k==1);

        [k1,thetaRefMajor,forceMajor, ...
         exMajor,eThetaMajor] = ...
            stageRHS( ...
                tk,Xk,K,cfg, ...
                previousMajorEx, ...
                previousMajorETheta, ...
                dt,initialMajor);

        thetaRef(k) = thetaRefMajor;
        force(k) = forceMajor;

        X2 = Xk+0.5*dt*k1;

        [k2,~,~,~,~] = ...
            stageRHS( ...
                tk+0.5*dt,X2,K,cfg, ...
                exMajor,eThetaMajor, ...
                0.5*dt,false);

        X3 = Xk+0.5*dt*k2;

        [k3,~,~,~,~] = ...
            stageRHS( ...
                tk+0.5*dt,X3,K,cfg, ...
                exMajor,eThetaMajor, ...
                0.5*dt,false);

        X4 = Xk+dt*k3;

        [k4,~,~,~,~] = ...
            stageRHS( ...
                tk+dt,X4,K,cfg, ...
                exMajor,eThetaMajor, ...
                dt,false);

        Xnext = ...
            Xk+dt*(k1+2*k2+2*k3+k4)/6;

        if any(~isfinite(Xnext)) || ...
           abs(Xnext(1))>cfg.maxAbsX || ...
           abs(Xnext(2))>cfg.maxAbsXd || ...
           abs(Xnext(3))>cfg.maxAbsTheta || ...
           abs(Xnext(4))>cfg.maxAbsThetaDot || ...
           abs(Xnext(5))>cfg.maxAbsIntegral || ...
           abs(Xnext(6))>cfg.maxAbsIntegral

            valid = false;
            break;
        end

        X(k+1,:) = Xnext.';

        previousMajorEx = exMajor;
        previousMajorETheta = eThetaMajor;
    end

    if ~valid

        sim = invalidSimulation(cfg,horizon);
        return;
    end

    [~,thetaRef(end),force(end),~,~] = ...
        stageRHS( ...
            t(end),X(end,:).',K,cfg, ...
            previousMajorEx, ...
            previousMajorETheta, ...
            dt,false);

    sim.t = t;
    sim.x = X(:,1);
    sim.xd = X(:,2);
    sim.theta = X(:,3);
    sim.thetad = X(:,4);
    sim.Ix = X(:,5);
    sim.Itheta = X(:,6);
    sim.thetaRef = thetaRef;
    sim.F = force;
end

function [dX,thetaRef,F,ex,eTheta] = ...
    stageRHS( ...
        t,X,K,cfg, ...
        derivativeBaseEx, ...
        derivativeBaseETheta, ...
        derivativeDeltaT, ...
        initialMajor)

    Kpx = K(1);
    Kix = K(2);
    Kdx = K(3);

    Kpp = K(4);
    Kip = K(5);
    Kdp = K(6);

    x = X(1);
    xd = X(2);
    theta = X(3);
    thetaDot = X(4);
    Ix = X(5);
    Itheta = X(6);

    %#ok<NASGU>
    t = t;

    ex = cfg.xRef-x;

    if initialMajor
        dEx = 0;
    else
        dEx = ...
            (ex-derivativeBaseEx) ...
            / derivativeDeltaT;
    end

    outerRaw = ...
        Kpx*ex ...
      + Kix*Ix ...
      + Kdx*dEx;

    outerLimited = clamp( ...
        outerRaw, ...
        -cfg.thetaRefSat, ...
        cfg.thetaRefSat);

    thetaRef = -outerLimited;

    eTheta = thetaRef-theta;

    if initialMajor
        dETheta = 0;
    else
        dETheta = ...
            (eTheta-derivativeBaseETheta) ...
            / derivativeDeltaT;
    end

    forceRaw = ...
        Kpp*eTheta ...
      + Kip*Itheta ...
      + Kdp*dETheta;

    F = clamp( ...
        forceRaw, ...
        -cfg.forceSat, ...
        cfg.forceSat);

    q = cfg.q;

    xdd = ...
        -(cfg.I+cfg.m*cfg.l^2)*cfg.b/q*xd ...
        + (cfg.m^2*cfg.g*cfg.l^2)/q*theta ...
        + (cfg.I+cfg.m*cfg.l^2)/q*F;

    thetaDD = ...
        -cfg.m*cfg.l*cfg.b/q*xd ...
        + cfg.m*cfg.g*cfg.l*(cfg.M+cfg.m)/q*theta ...
        + cfg.m*cfg.l/q*F;

    IxDot = ex;
    IthetaDot = eTheta;

    dX = [ ...
        xd;
        xdd;
        thetaDot;
        thetaDD;
        IxDot;
        IthetaDot];
end

function [met,ok] = ...
    calculateMetrics(sim,cfg,horizon)

    ok = ...
        all(isfinite(sim.x)) && ...
        all(isfinite(sim.theta)) && ...
        all(isfinite(sim.F));

    if ~ok
        met = invalidMetrics();
        return;
    end

    met.thetaMax = max(abs(sim.theta));
    met.thetaRefMax = max(abs(sim.thetaRef));
    met.Fmax = max(abs(sim.F));

    met.forceSatRatio = ...
        mean(abs(sim.F)>=0.99*cfg.forceSat);

    met.thetaRefSatRatio = ...
        mean(abs(sim.thetaRef)>=0.99*cfg.thetaRefSat);

    tail = sim.t>=0.8*horizon;

    met.xTailStd = std(sim.x(tail));
    met.thetaTailStd = std(sim.theta(tail));

    tailSS = sim.t>=0.9*horizon;

    met.essX = abs( ...
        cfg.xRef-mean(sim.x(tailSS)));

    met.OSx = max( ...
        0,max(sim.x)-cfg.xRef);

    met.OSxPercent = ...
        100*met.OSx/max(abs(cfg.xRef),eps);

    met.TsX = settlingTime( ...
        sim.t,sim.x,cfg.xRef, ...
        cfg.xTolerance);

    met.TsTheta = settlingTime( ...
        sim.t,sim.theta,0, ...
        cfg.thetaTolerance);

    met.TsOverall = max( ...
        met.TsX,met.TsTheta);

    met.pass = ...
        isfinite(met.TsOverall) && ...
        met.TsOverall<=cfg.TsReq && ...
        met.essX<=cfg.xEssReq && ...
        met.OSx<=cfg.maxOvershoot && ...
        met.thetaMax<=cfg.thetaLimit && ...
        met.Fmax<=cfg.forceSat+1e-9;

    met.violation = inf;
    met.quality = inf;
    met.maxPoleReal = inf;
    met.poleMarginPass = false;
    met.poleMarginViolation = inf;
    met.poles = NaN(6,1);
    met.characteristicPolynomial = NaN(1,7);
end

function Ts = settlingTime( ...
    t,y,reference,tolerance)

    lastOutside = find( ...
        abs(y-reference)>tolerance, ...
        1,'last');

    if isempty(lastOutside)

        Ts = 0;
        return;
    end

    if lastOutside>=numel(t)

        Ts = inf;
    else
        Ts = t(lastOutside+1);
    end
end

function result = emptyResult(cfg,horizon)

    result.K = ones(1,6);
    result.fitness = inf;
    result.met = invalidMetrics();
    result.sim = invalidSimulation(cfg,horizon);
end

function printResult(result)

    fprintf('Fitness     = %.10g\n', ...
        result.fitness);
    fprintf('Feasible    = %d\n', ...
        result.met.pass);
    fprintf('Violation   = %.10g\n', ...
        result.met.violation);
    fprintf('Ts overall  = %.8g s\n', ...
        result.met.TsOverall);
    fprintf('Overshoot   = %.8g %%\n', ...
        result.met.OSxPercent);
    fprintf('ess         = %.8g m\n', ...
        result.met.essX);
    fprintf('theta max   = %.8g deg\n', ...
        rad2deg(result.met.thetaMax));
    fprintf('F max       = %.8g N\n', ...
        result.met.Fmax);
    fprintf('Max Re pole = %.8g 1/s\n', ...
        result.met.maxPoleReal);
    fprintf('Pole robust = %d\n', ...
        result.met.poleMarginPass);
end

function z = guardLogVector(z,maxAbsValue)

    z = min(max(z,-maxAbsValue),maxAbsValue);
end

function y = clamp(u,lowerLimit,upperLimit)

    y = min(max(u,lowerLimit),upperLimit);
end

function sim = invalidSimulation(cfg,horizon)

    t = (0:cfg.dt:horizon).';
    n = numel(t);

    sim.t = t;
    sim.x = nan(n,1);
    sim.xd = nan(n,1);
    sim.theta = nan(n,1);
    sim.thetad = nan(n,1);
    sim.Ix = nan(n,1);
    sim.Itheta = nan(n,1);
    sim.thetaRef = nan(n,1);
    sim.F = nan(n,1);
end

function met = invalidMetrics()

    met.TsX = inf;
    met.TsTheta = inf;
    met.TsOverall = inf;

    met.essX = inf;
    met.OSx = inf;
    met.OSxPercent = inf;

    met.thetaMax = inf;
    met.thetaRefMax = inf;
    met.Fmax = inf;

    met.forceSatRatio = 1;
    met.thetaRefSatRatio = 1;

    met.xTailStd = inf;
    met.thetaTailStd = inf;

    met.pass = false;

    met.violation = inf;
    met.quality = inf;
    met.maxPoleReal = inf;
    met.poles = NaN(6,1);
    met.characteristicPolynomial = NaN(1,7);
end
