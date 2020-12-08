% Runs a Genetic Algorithm on the objective function f(x) with constraint
% space [[0;0;0], [20;20;20]].
function GAOptimize(P, S, ITL, enc_prec, gamma)
	% Seeding the random number generator to ensure it generates the same
    % random numbers each time.
    rng(0, 'twister');
    % Used for the file output.
    % The algorithm aims to record 25 data points
    logi = max(floor(ITL / 25), 1);
    % Log File
    file_fbest = fopen("ga_fbest.csv", "w");    
    % Setting up n and the contraint set [l, u]
    n = 3;
    l = [0; 0; 0];
    u = [20; 20; 20];
    % Allocating space for the P points in each population
    Population = zeros(n, P);
    % d represents the number of bits used to encode each coordinate.
    % It is calculated using the encoding precision variable. The smaller
    % enc_prec is, the larger d will be.
    d = 0;
    for i = 1 : n
        d = max(d, ceil(log2((u(i) - l(i))/enc_prec)));
        % Generating starting population
        Population(i, :) = l(i) + (u(i) - l(i))*rand(1, P);
    end
    % Variables to keep track of the best function value found and what 
    % input produces this function value.
    xbest = Population(:, 1);
    fbest = Inf;
    % Interval between encoded points
    Delta = (u - l) / pow2(d);
    % This loop runs ITL times. Each time it runs, it evaluates f at each
    % point in a population. Using this information, a new population is
    % generated.
    for p = 1 : ITL
        % Evaluating f at each point in the population
        f_evals = zeros(1, P);
        for y = 1 : P
            f_evals(y) = f(Population(:, y));
            % Updating best function value and input
            if f_evals(y) < fbest
                xbest = Population(:, y);
                fbest = f_evals(y);
            end
        end
        % Output to the user
        fprintf("Population %d\n\t- f_best = %g\n\n", p, fbest)
        % Output to log every logi iterations
        if mod(p, logi) == 0
            fprintf(file_fbest, "%g, ", fbest);
        end
        % Computes the fitness of each population member. This can be
        %   - fit_rank(f_evals, P)
        %   - fit_value(f_evals, P)
        fitness = fit_rank(f_evals, P);
        % Allocating space for the next population
        nPop = zeros(n, P);
        % Copying the fitness list to use for selecting survivors. So:
        %   - tmp_fit will be used for selecting survivors
        %   - fitness will be used for selecting parents
        tmp_fit = fitness;
        % Filling the new population list with S survivors
        for s = 1 : S
            % Selects the next survivor. This can be
            %   - selection_elitism(tmp_fit, P)
            %   - selection_roulette(tmp_fit, P)
            %   - selection_etournament(tmp_fit, P, m)
            %   - selection_rtournament(tmp_fit, P, m)
            % where m indicates the tournament size.
            surv = selection_roulette(tmp_fit, P);
            % Copying the selected survivor to the new population list
            nPop(:, s) = Population(:, surv);
            % This ensures that no survivor is picked more than once.
            tmp_fit(surv) = 0;
        end
        
        % Generating offspring
        % This loop runs until nPop is filled completely
        o = S + 1;
        while o <= P
            % Selecting both parents (similar to survivor selection).
            % However, fitness is used instead of tmp_fit so that survivors
            % can also be parents.
            parent1 = selection_roulette(fitness, P);
            fitness(parent1) = 0;
            parent2 = selection_roulette(fitness, P);
            fitness(parent2) = 0;
            % Encoding each parent using binary discretisation
            gen_p1 = bin_enc(Population(:, parent1), l, d, Delta, n);
            gen_p2 = bin_enc(Population(:, parent2), l, d, Delta, n);
            % Generating the child genome. This can be
            %   - crossover_1point(gen_p1, gen_p2, n)
            %   - crossover_2point(gen_p1, gen_p2, n)
            %   - crossover_probselect(gen_p1, gen_p2, n, m)
            % where m = P(gen_p1(i) = gen_child(i)).
            gen_child = crossover_1point(gen_p1, gen_p2, n);
            % Mutates the child genome with probability gamma
            if rand() < gamma
                gen_child = mutation_binv(gen_child, n, d, 0.1);
            end
            % Convert the child genome into the actual coordinates.
            child = bin_dec(gen_child, l, Delta, n);
            % Check if the coordinates are still in the constraint set
            alive = true;
            for i = 1 : n
                if ~(l(i) <= child(i) <= u(i))
                    alive = false;
                end
            end
            % If the child is valid, it is added to the new population list
            if alive
                nPop(:, o) = child;
                o = o + 1;
            end
        end
        % After the new population is fully generated, it replaces the old
        % population
        Population = nPop;
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

% Returns a list containing the rank fitness of each f-value in f_evals.
function fitness = fit_rank(f_evals, P)
    fitness = zeros(1, P);
    for y = 1 : P
        for z = 1 : P
            if f_evals(y) <= f_evals(z)
                fitness(y) = fitness(y) + 1;
            end
        end
    end
end

% Returns a list containing the function value fitness of each f-value in
% f_evals.
function fitness = fit_value(f_evals, P)
    fitness = zeros(1, P);
    f_bar = max(f_evals);
    for y = 1 : P 
        fitness(y) = f_bar - f_evals(y) + 1;
    end
end

% Returns the index of the population member with the highest fitness.
function S = selection_elitism(fitness, P)
    f_bar = max(fitness);
    for i = 1 : P
        if fitness(i) == f_bar
            S = i;
            break;
        end
    end
end

% Returns a population member that is picked using roulette selection.
function S = selection_roulette(fitness, P)
    % Calculating the probability for each population member.
    f_sum = sum(fitness);
    probs = zeros(1, P);
    for i = 1 : P
        probs(i) = (fitness(i) / f_sum);
        if i > 1
            probs(i) = probs(i) + probs(i - 1);
        end
    end
    % Picking the winner
    roll = rand(1);
    winner = 0;
    for i = 1 : P
        if roll < probs(i)
            winner = i;
            break;
        end
    end
    S = winner;
end

% Returns a population member that is picked using tournament selection.
% The tournament is of size m and is decided using elitism selection.
function S = selection_etournament(fitness, P, m)
    % Generating a tournament of size m (this may be less if less than m
    % member remain in fitness).
    nz_elem = size(nonzeros(fitness));
    m = min(m, nz_elem(1));
    tournaments = zeros(2, m);
    for i = 1 : m
        roll = uint8(rand(1) * (P - 1) + 1);
        while ismember(roll, tournaments(1, :)) || fitness(roll) == 0
            roll = mod(roll, P) + 1;
        end
        tournaments(1, i) = roll;
        tournaments(2, i) = fitness(roll);
    end
    % Using elitism selection on the generated tournament
    t_winner = selection_elitism(tournaments(2, :), m);
    S = tournaments(1, t_winner);
end

% Returns a population member that is picked using tournament selection.
% The tournament is of size m and is decided using roulette selection.
function S = selection_rtournament(fitness, P, m)
    % Generating a tournament of size m (this may be less if less than m
    % member remain in fitness).
    nz_elem = size(nonzeros(fitness));
    m = min(m, nz_elem(1));
    tournaments = zeros(2, m);
    for i = 1 : m
        roll = uint8(rand(1) * (P - 1) + 1);
        while ismember(roll, tournaments(1, :)) || fitness(roll) == 0
            roll = mod(roll, P) + 1;
        end
        tournaments(1, i) = roll;
        tournaments(2, i) = fitness(roll);
    end
    % Using roulette selection on the generated tournament
    t_winner = selection_roulette(tournaments(2, :), m);
    S = tournaments(1, t_winner);
end

% Encodes a point in the constraint space using binary discretisation.
function xenc = bin_enc(x, l, d, Delta, n)
    % For each coordinate, xi in [li + mi Deltai, li + (mi + 1)Deltai), the
    % value mi is calculated and stored in the list m.
    m = zeros(n, 1);
    for i = 1 : n
        m(i) = floor((x(i) - l(i))/Delta(i));
        if m(i) == pow2(d)
            m(i) = m(i) - 1;
        end
    end
    % The list m is then converted into a matrix where each row is mi in
    % binary.
    xenc = de2bi(m, d, 'left-msb');
end

% Decodes a point in the constraint space using binary discretisation.
function xdec = bin_dec(x, l, Delta, n)
    % The encoding is reversed to restore the list m (see bin_enc).
    k = bi2de(x, 'left-msb');
    % To generate a point in Rn, the midpoint of each interval is used.
    xdec = zeros(n, 1);
    for i = 1 : n
        xdec(i) = l(i) + (k(i) + 0.5)*Delta(i);
    end
end

% Generates a child using single-point crossover and with the parent
% genomes p1, p2.
function child = crossover_1point(p1, p2, n)
    child = p1;
    X = randi([1 (n - 1)]);
    
    child((X+1):n, :) = p2((X+1):n, :);
end

% Generates a child using 2-point crossover and with the parent genomes p1, 
% p2.
function child = crossover_2point(p1, p2, n)
    child = p1;
    X = randi([1 (n-2)]);
    Y = randi([(X + 1) (n - 1)]);
    
    child((X+1):Y, :) = p2((X+1):Y, :);
end

% Generates a child using probability select crossover and with the parent
% genomes p1, p2.
function child = crossover_probselect(p1, p2, n, theta)
    child = p1;
    for r = 1 : n
        if rand() > theta
            child(r, :) = p2(r, :);
        end
    end
end

% Generates a new genome given gen using bit inversion.
function child = mutation_binv(gen, n, d, delta)
    child = gen;
    for r = 1 : n
        for c = 1 : d
            if rand() < delta
                if child(r, c) == 1
                    child(r, c) = 0;
                else
                    child(r,c) = 1;
                end
            end
        end
    end
end