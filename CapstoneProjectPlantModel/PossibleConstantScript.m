clear;
clc;
warning('off');
Cascade_SolutionPumpModel_DD;

PumpPressureConstant = 3.275E-5; % Calculated from previous work in Simulink

% Load the Simulink model
load_system('BasicPumpPlant.slx');

if ~isfile("pump_constants.mat")

    % Desired range of target pressures and target flows
    target_pressure_range = 10:10:400; % psi (example: 10 to 400 psi, in steps of 10)
    target_flow_range = 0:10:800; % gpm (example: 0 to 800 gpm, in steps of 10)

    % Initialize matrix to store results (rows: target pressures, columns: target flows)
    speed_results = zeros(length(target_pressure_range), length(target_flow_range));
    pressure_results = zeros(length(target_pressure_range), length(target_flow_range));
    flow_results = zeros(length(target_pressure_range), length(target_flow_range));
    dead_head_pressure = zeros(length(target_pressure_range), length(target_flow_range));
    PumpCoeff = zeros(length(flow_results), length(flow_results));

    % Loop through target pressures and target flows
    total_entries = (length(target_pressure_range)) * (length(target_flow_range));
    f = waitbar(0, "Simulating...");
    for i = 1:length(target_pressure_range)
        for j = 1:length(target_flow_range)
            % Print a progress bar so we know how far we are from
            % completing data collection
            percent_complete = ((i-1)*length(target_flow_range) + j)/total_entries;
            waitbar(percent_complete, f, "Simulating: " + round(percent_complete*100) + "%")

            % Set target pressure and target flow in the BasicPumpPlant model
            set_param('BasicPumpPlant/TP_Constant', 'Value', num2str(target_pressure_range(i)));
            set_param('BasicPumpPlant/TF_Constant', 'Value', num2str(target_flow_range(j)));

            % Run the Simulink model
            sim('BasicPumpPlant.slx');

            % Record the resulting speed, pressure, and flow from the SolutionPumpPlantModel
            speed_results(i, j) = ShaftSpeedModel_rpm.signals.values(end);
            pressure_results(i, j) = PumpPressureModel_psi.signals.values(end);
            flow_results(i, j) = PumpFlowModel_gpm.signals.values(end);
%         end
%     end
%     waitbar(1,f,'Finishing');
%     pause(1)
%     close(f)


% Find the dead head pressure - Scalar (I think this equation is wrong but I'm not sure what it should be)
% Dead head pressure is with target flow == 0, So should this 
        dead_head_pressure(i, j) = 6.89476 * PumpPressureConstant / (pressure_results(i,1)^2);

        end
    end
    waitbar(1,f,'Finishing');
    pause(1)
    close(f)
else
    sim('BasicPumpPlant.slx');
    load("pump_constants.mat");
end
% Calculate the pressure drop (matrix values at each target pressure/flow)
pressure_drop = dead_head_pressure - pressure_results;

% Fit a polynomial to the pressure drop data
% Order 10 is the only thing that gets a non-negative PumpACd (Red Flag)
poly_coeff = polyfit(flow_results, pressure_drop, 10); % Quadratic polynomial fit

% Calculate the PumpACdConstant
PumpACdConstant = 1 / sqrt(poly_coeff(1));

% Create target pressure and target flow matrices
[TargetPressureMesh, TargetFlowMesh] = meshgrid(target_pressure_range, target_flow_range);

% Reshape the matrices to create lists for plotting
TargetPressureVec = reshape(TargetPressureMesh, [], 1);
TargetFlowVec = reshape(TargetFlowMesh, [], 1);
PressureResultsVec = reshape(pressure_results', [], 1);

% Create a 3D scatter plot
figure
scatter3(TargetPressureVec, TargetFlowVec, PressureResultsVec, 'filled')
xlabel('Target Pressure (psi)')
ylabel('Target Flow (gpm)')
zlabel('Pressure Results (psi)')
title('3D Scatter Plot of Pressure Results')
grid on

% Save the results and constants to a .mat file
save('pump_constants.mat', 'speed_results', 'pressure_results', 'flow_results', 'target_pressure_range', 'target_flow_range', 'PumpACdConstant');