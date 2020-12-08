% Runs the GS Algorithm on the interval [[0;0;0], [20;20;20]] to minimize
% f(x).
function GSOptimize(p)
    % Used for the file output.
    % The algorithm aims to record 25 data points
    logi = max(floor((p^3) / 25), 1);
    % Setting up n and the contraint set [l, u]
    n = 3;
    l = [0;0;0];
    u = [20;20;20];
    % Log File
    file_fbest = fopen("gs_fbest.csv", "w");
    % Size of each interval between points
    delta = (u - l) / (p - 1);
    % Variables to keep track of the best function value found and what 
    % input produces this function value.
    xbest = l;
    fbest = Inf;
    % Iteration counter for the purposes of logging fbest
    IT = 1;
    % Trying every point in the grid as an input value for f(x).
    for a = 0 : p - 1
        for b = 0 : p - 1
            for c = 0 : p - 1
                % Constructing xk and finding fk = f(xk)
                xk = zeros(n, 1);
                xk(1) = l(1) + a * delta(1);
                xk(2) = l(2) + b * delta(2);
                xk(3) = l(3) + c * delta(3);
                fk = f(xk);
                % Updating best function value and input
                if fk < fbest
                    xbest = xk;
                    fbest = fk;
                end
                % Output to log every logi iterations
                if mod(IT, logi) == 0
                    fprintf(file_fbest, "%g, ", fbest);
                end
                IT = IT + 1;
            end
        end
    end
    % Closing the log file
    fclose(file_fbest);
    % Displaying result to user
    fprintf("\n\nAlgorithm Finished\n\nf_best = %g\n", fbest)
    fprintf("with x_best = [%g %g %g]^T\n", xbest(1), xbest(2), xbest(3));
end

% This version squares the error term (ie: z = ∑ε(x)^2)
function z = g_smooth(data, x)
    z = 0;
    for i = 1 : 13
        z = z + (abs(x(1) * (1 + (x(2)^2) * (data(i, 1)^2))^((x(3) - 1) / 2) - data(i, 2)))^2;
    end
end

% This version just adds up each error value (ie: z = ∑ε(x))
function z = g_nonsmooth(data, x)
    z = 0;
    for i = 1 : 13
        z = z + (abs(x(1) * (1 + (x(2)^2) * (data(i, 1)^2))^((x(3) - 1) / 2) - data(i, 2)));
    end
end

% Objective Function
function z = f(x)
    data = [0.0137, 3220;0.0274, 2190;0.0434, 1640;0.0866, 1050;0.137, 766;0.274, 490;0.434, 348;0.866, 223;1.37, 163;2.74, 104;4.34, 76.7;5.46, 68.1;6.88, 58.2];
    z = g_nonsmooth(data, [520 * x(1); 14 * x(2); 0.038 * x(3)]);
end