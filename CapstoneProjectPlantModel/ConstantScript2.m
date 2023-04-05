clear;
clc;
% warning('off');
Cascade_SolutionPumpModel_DD;

PumpPressureConstant = 3.275E-5; % Calculated from previous work in Simulink

% Load the Simulink model
load_system('BasicPumpPlant.slx');

if ~isfile("pump_constants.mat")

    % Desired range of target pressures and target flows
    input_speed_current = 10:500:6000; % psi (example: 10 to 400 psi, in steps of 10)
    input_DSOutletArea = 0:.5:5; % m^2 (example: 0 to 800 m^2, in steps of 10) 2.5" plumbing - Area is 4.91

    % Initialize matrix to store results (rows: target pressures, columns: target flows)
    speed_results = zeros(length(input_speed_current), length(input_DSOutletArea));
    pressure_results = zeros(length(input_speed_current), length(input_DSOutletArea));
    flow_results = zeros(length(input_speed_current), length(input_DSOutletArea));
    dead_head_pressure = zeros(length(input_speed_current), length(input_DSOutletArea));
    PumpCoeff = zeros(length(flow_results), length(flow_results));

    % Loop through target pressures and target flows
    % ** Loop through speed and dsoutletarea to characterize the pump
    total_entries = (length(input_speed_current)) * (length(input_DSOutletArea));
    f = waitbar(0, "Simulating...");
    for i = 1:length(input_speed_current)
        for j = 1:length(input_DSOutletArea)
            % Print a progress bar so we know how far we are from
            % completing data collection
            percent_complete = ((i-1)*length(input_DSOutletArea) + j)/total_entries;
            waitbar(percent_complete, f, "Simulating: " + round(percent_complete*100) + "%  Speed: " + input_speed_current(i) + " DSOA: " + input_DSOutletArea(j))

            % Set target pressure and target flow in the BasicPumpPlant model
            set_param('BasicPumpPlant/TS_Constant', 'Value', num2str(input_speed_current(i)));
            set_param('BasicPumpPlant/TDSOA_Constant', 'Value', num2str(input_DSOutletArea(j)));

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
% **** Should these be flipped? ****
pressure_drop = dead_head_pressure - pressure_results;

% Fit a polynomial to the pressure drop data
% Order 10 is the only thing that gets a non-negative PumpACd (Red Flag)
poly_coeff = polyfit(flow_results, pressure_drop, 2); % Quadratic polynomial fit
Pressure_coeff = poly_coeff(3) ./ (speed_results.^2); % Use 

% Calculate the PumpACdConstant
%  PumpACdConstant = 1 / sqrt(poly_coeff(2));
PumpACdConstant = 1 ./ sqrt(Pressure_coeff);
display(PumpACdConstant);

% Use Peterman equation to solve for Impeller Drag Algebraicly

% Create target pressure and target flow matrices
[InputSpeedMesh, InputDSOAMesh] = meshgrid(input_speed_current, input_DSOutletArea);
% 
% % Reshape the matrices to create lists for plotting
InputSpeedVec = reshape(InputSpeedMesh, [], 1);
InputDSOAVec = reshape(InputDSOAMesh, [], 1);
PressureResultsVec = reshape(pressure_results', [], 1);
FlowResultsVec = reshape(flow_results', [], 1);

% Create a 3D scatter plot
figure
scatter3(InputSpeedVec, InputDSOAVec, PressureResultsVec, 'filled')
xlabel('Input Speed')
ylabel('Input DSOA')
zlabel('Pressure Results (psi)')
title('3D Scatter Plot of Pressure Results')
grid on

% Create a 3D scatter plot
figure
scatter3(InputSpeedVec, InputDSOAVec, FlowResultsVec, 'filled')
xlabel('Input Speed')
ylabel('Input DSOA')
zlabel('Flow Results')
title('3D Scatter Plot of Pressure Results')
grid on

figure
scatter(FlowResultsVec, PressureResultsVec, 'filled')
ylabel('Pressure Results (psi)')
xlabel('Flow Results (gpm)')
title('2D Scatter Plot of Pump Curve')
grid on

% Save the results and constants to a .mat file
save('pump_constants.mat', 'speed_results', 'pressure_results', 'flow_results', 'input_speed_current', 'input_DSOutletArea', 'dead_head_pressure', 'PumpACdConstant');