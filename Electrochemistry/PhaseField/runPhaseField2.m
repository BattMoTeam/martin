%% Read model parameters
%

filename = fullfile(battmoDir(), 'Electrochemistry', 'PhaseField', 'jsonfiles', 'phasefield.json');
jsonstruct = parseBattmoJson(filename);
jsonstruct.np = 3;
    
% instantiate model (done also by setupPhaseFieldSimulation)
inputparams = PhaseFieldInputParams(jsonstruct);
model = PhaseField(inputparams);

%% Setup model, initial state

simsetup = setupPhaseFieldSimulation(jsonstruct);

%% Setup schedule
% in this case only the time steps are given, no source term

total = 5;
n     = 100;                       
dt    = total / n;                  
dts   = rampupTimesteps(total, dt, 5);

% define flux boundary condition function during the time of the simulation
% the sign convension is that the outward flux is positive.
clear flux
flux.functionFormat = 'tabulated';
flux.argumentList   = {'time'};
% flux.dataX          = [0, 10, 11, 20, 21, total];
% flux.dataY          = [0, 0, 0.01, 0.01, -0.01, -0.01];
flux.dataX = [0, total];
flux.dataY = [0, 0];

fluxFunc = setupFunction(flux);
control  = struct('src', @(t) fluxFunc(t)*ones(model.np, 1));

step = struct('val', dts, 'control', ones(numel(dts), 1));
schedule = struct('control', control, 'step', step);

simsetup.schedule = schedule;

%% Run simulation
simsetup.model.nonlinearTolerance = 1e-20;
simsetup.model.verbose = true;
% simsetup.run();


%% Visualization
[states, globvars, reports] = simsetup.run();
model = simsetup.model;

%%

x = model.chebyshevNodes(model.N);  
[c1, c2] = model.getEquilibriumValues(1/PhysicalConstants.kb);


% %% Figure 1 : energy functions
% c_test = linspace(0.01, 0.99, 200)';
% 
% % F'(c) used in eqW
% f_prime = model.energyFunc(c_test);
% 
% % f_hom(c) double well shape
% omega = model.omega;
% kT    = model.kT;
% f_hom = omega .* c_test .* (1 - c_test) + kT .* (c_test .* log(c_test) + (1 - c_test) .* log(1 - c_test));
% 
% figure;
% 
% subplot(1, 2, 1);
% plot(c_test, f_prime, 'b-', 'LineWidth', 1.5);
% xlabel('c');
% ylabel("F'(c)");
% title('Energy derivative (used in eqW)');
% grid on;
% 
% subplot(1, 2, 2);
% plot(c_test, f_hom, 'r-', 'LineWidth', 1.5);
% xlabel('c');
% ylabel('f_{hom}(c)');
% title('Free energy : double well');
% grid on;



%% Figure 1.5 : boundary flux control over time
t_plot     = linspace(0, total, 500);
fluxValues = arrayfun(@(t) fluxFunc(t), t_plot);

figure;
plot(t_plot, fluxValues, 'b-', 'LineWidth', 1.5);
xlabel('time');
ylabel('J_{in}');
title('Boundary flux control at x = 1');
grid on;
yline(0, 'k--', 'LineWidth', 0.5);


% %% Figure 2 : concentration profiles
% nPlots = numel(dts);
% nPlots = 10;
% nPlots = min(nPlots, numel(states));  % in case fewer states are saved
% 
% figure;
% hold on;
% 
% % initial state with black line
% cInit = simsetup.initstate.c;
% plot(x, cInit, 'k-', 'LineWidth', 2, 'DisplayName', 'initial');
% 
% % states 1 to nPlots with color gradient
% cmap = nebula(nPlots);
% for iState = 1 : nPlots
%     c = states{iState}.c;
%     plot(x, c, 'Color', cmap(iState, :), 'DisplayName', sprintf('state %d', iState));
% end
% 
% if nPlots <= 10
%     legend('show', 'Location', 'best');
% end
% colormap('nebula')
% cb = colorbar();
% cb.Ticks      = [0, 1];
% cb.TickLabels = {'Initial state', 'Final state'};
% xlabel('x');
% ylabel('c');
% title(sprintf('Concentration profiles : initial + states 1 to %d', nPlots));
% grid on;
% hold off;




% %% Animated plot of concentration profile
% times = cellfun(@(s) s.time, states);
% 
% figure;
% hLine = plot(x, model.A * simsetup.initstate.coefC, 'k-', 'LineWidth', 2);
% xlabel('x');
% ylabel('c');
% title('Concentration profile c(x,t)');
% grid on;
% ylim([0, 1]);  % fix y axis to avoid rescaling at each frame
% 
% nStates = numel(states);
% nTrail  = 20;  % number of trailing curves to keep
% cmap = turbo(nTrail);
% 
% for iState = 1 : nStates
%     c = states{iState}.c;
% 
%     % update the existing line instead of redrawing
%     set(hLine, 'YData', c);
%     title(sprintf('Concentration profile c(x,t)  —  t = %.4f  (state %d / %d)', ...
%                   times(iState), iState, nStates));
% 
%     drawnow;
%     pause(0.005);  % adjust to control animation speed
% end


% 
% %% Animation with trail and color gradient
% 
% times   = cellfun(@(s) s.time, states);
% nStates = numel(states);
% nTrail  = 15;  % number of trailing curves to keep
% 
% figure;
% ax = axes;
% xlabel(ax, 'x');
% ylabel(ax, 'c');
% grid(ax, 'on');
% ylim(ax, [0, 1]);
% hold(ax, 'on');
% 
% cmap = turbo(nTrail);
% for iState = 1 : nStates
% 
%     cla(ax);  % clear axes but keep settings
% 
%     % draw trail : last nTrail states
%     iStart = max(1, iState - nTrail + 1);
%     nDrawn = iState - iStart + 1;
% 
%     for iTrail = iStart : iState
%         alpha     = (iTrail - iStart + 1) / nDrawn;  % 0 = oldest, 1 = newest
% 
%         trailIdx  = round(alpha * (nTrail - 1)) + 1;
%         c         = states{iTrail}.c;
%         lineWidth = 0.5 + 1.5 * alpha;  % thicker for more recent
%         plot(ax, x, c, 'Color', [cmap(trailIdx, :), alpha], 'LineWidth', lineWidth);
%     end
% 
%     title(ax, sprintf('c(x,t)  —  t = %.4f  (state %d / %d)', ...
%                       times(iState), iState, nStates));
%     drawnow;
%     pause(0.02);
% end


% %% Animation with color based on control phase
% % phase 1 : t in [0, 10]   -> no flux     -> blue
% % phase 2 : t in [11, 20]  -> flux = 0.01 -> red
% % phase 3 : t in [21, end] -> flux = -0.01 -> green
% 
% function col = getPhaseColor(t)
%     if t <= 10
%         col = [0.2, 0.4, 0.8];   % blue  — no flux
%     elseif t <= 20
%         col = [0.8, 0.2, 0.2];   % red   — positive flux
%     else
%         col = [0.2, 0.7, 0.3];   % green — negative flux
%     end
% end
% 
% figure;
% ax = axes;
% xlabel(ax, 'x');
% ylabel(ax, 'c');
% grid(ax, 'on');
% ylim(ax, [0, 1]);
% hold(ax, 'on');
% 
% for iState = 1 : nStates
%     cla(ax);
%     iStart = max(1, iState - nTrail + 1);
%     for iTrail = iStart : iState
%         alpha    = (iTrail - iStart + 1) / (iState - iStart + 1);
%         col      = getPhaseColor(times(iTrail));
%         c        = states{iTrail}.c;
%         plot(ax, x, c, 'Color', [col, alpha], 'LineWidth', 0.5 + 1.5*alpha);
%     end
%     title(ax, sprintf('c(x,t)  —  t = %.4f  (state %d / %d)', ...
%                       times(iState), iState, nStates));
%     drawnow;
%     pause(0.003);
% end



%% Interactive visualization with slider for multiples particles

np = model.np;
times   = cellfun(@(s) s.time, states);
nStates = numel(states);

% vis parameters into a struct
op = model.operators;

plotParams.x            = x;
plotParams.nTrail       = 1;  % number of trailing curves
plotParams.cmap         = turbo(plotParams.nTrail);
plotParams.c0           = simsetup.initstate.c;
plotParams.indInnerBc   = op.indInnerBc;
plotParams.nPt          = model.N + 1;
plotParams.np           = np;
plotParams.ipSelected = 1;  % default : show all particles

% create figure 
fig = figure;
ax  = axes(fig, 'Position', [0.1, 0.2, 0.85, 0.75]);
xlabel(ax, 'x');
ylabel(ax, 'c');
grid(ax, 'on');
ylim(ax, [0, 1]);
hold(ax, 'on');


% time slider
sld = uicontrol(fig, 'Style', 'slider', ...
                'Min', 1, 'Max', nStates, 'Value', 1, ...
                'SliderStep', [1/(nStates-1), 10/(nStates-1)], ...
                'Position', [80, 20, 640, 20]);

% time label 
lbl = uicontrol(fig, 'Style', 'text', ...
                'Position', [80, 45, 640, 20], ...
                'String', 'state 1');

% particle selection buttons
bgParticle = uibuttongroup(fig, ...
    'Title',    'Particle', ...
    'Position', [0.02, 0.15, 0.05, 0.75]);

% 'all' button at the top
uicontrol(bgParticle, 'Style', 'radiobutton', ...
          'String',   'all', ...
          'Tag',      '0', ...
          'Position', [5, np * 28 + 5, 70, 22], ...
          'Value',    1);

% one button per particle
for ip = 1 : np
    uicontrol(bgParticle, 'Style', 'radiobutton', ...
              'String',   sprintf('p %d', ip), ...
              'Tag',      sprintf('%d', ip), ...
              'Position', [5, (np - ip) * 28 + 5, 70, 22], ...
              'Value',    0);
end

% particle selection callback
bgParticle.SelectionChangedFcn = @(src, ~) updatePlot(sld, ax, states, times, lbl, plotParams, src);

% slider handler 
addlistener(sld, 'ContinuousValueChange', @(src, ~) updatePlot(src, ax, states, times, lbl, plotParams, bgParticle));


% draw first step
updatePlot(sld, ax, states, times, lbl, plotParams, bgParticle);

function updatePlot(sld, ax, states, times, lbl, plotParams, bgParticle)

    % check that all graphics objects are still valid
    if ~isvalid(sld) || ~isvalid(ax) || ~isvalid(lbl)
        return;
    end

    % retrieve plot params
    x            = plotParams.x;
    nTrail       = plotParams.nTrail;  
    cmap         = plotParams.cmap;
    c0           = plotParams.c0;
    indInnerBc   = plotParams.indInnerBc;
    nPt          = plotParams.nPt;
    np           = plotParams.np;

    % get selected particle from button group 
    ipTag      = str2double(bgParticle.SelectedObject.Tag);
    if ipTag == 0
        ipList = 1 : np;  % show all particles
    else
        ipList = ipTag;   % show only selected particle
    end

    
    iState = round(sld.Value);
    cla(ax); % clear axes but keep settings
    hold(ax, 'on');

    % keep initial state displayed
    for ip = ipList
        idx = indInnerBc(ip) : indInnerBc(ip) + nPt - 1;
        plot(ax, x, c0(idx), 'k-', 'LineWidth', 1, 'DisplayName','initial');
    end
    
    % main plot
    iStart = max(1, iState - nTrail + 1);
    nDrawn = iState - iStart + 1;
    for iTrail = iStart : iState
        alpha    = (iTrail - iStart + 1) / nDrawn; % 0 = oldest, 1 = newest
        trailIdx = round(alpha * (nTrail - 1)) + 1;
        c        = states{iTrail}.c;
        for ip = ipList
            idx = indInnerBc(ip) : indInnerBc(ip) + nPt - 1;
            plot(ax, x, c(idx), ...
                 'Color',     [cmap(trailIdx, :), alpha], ...
                 'LineWidth', 0.5 + 1.5 * alpha);
        end
    end

    title(ax, sprintf('c(x,t)  @  t = %.2f  (state %d / %d)', times(iState), iState, numel(states)));
    lbl.String = sprintf('state %d / %d  @  t = %.4f', iState, numel(states), times(iState));
    % ylim(ax, [0, 1]);
    grid(ax, 'on');

end
