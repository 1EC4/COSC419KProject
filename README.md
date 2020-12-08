## Running the Algorithms

1. Download one or more of the above `.m` files.
2. Ensure that MATLAB's working directory is set to the same directory containing the `.m` files you downloaded.
3. Execute the `.m` file.

#### GAOptimize

```
> GAOptimize(P, S, ITL, enc_prec, gamma)
```
where
- `P` is the size of each population.
- `S` is the number of survivors selected from each population.
- `ITL` is the total number of populations to generate.
- `enc_prec` represents the encoding precision. It is the maximum distance between points after encoding and decoding them.
- `gamma` is the probability that a child is mutated before being added to the a new population.


#### GSOptimize

```
> GSOptimize(p)
```
where `p` is the number of point *per dimension*. So for `n = 3`, there will be `p^n` points.


#### LHSOptimize

```
> LHSOptimize(p)
```
where `p` is the number of points to select using the LHS algorithm.


## Customizing GAOptimize

1. The function that generates the fitness list can be:
   - `fit_rank(f_evals, P);`
   - `fit_value(f_evals, P);`
   
   Replace the function on line `56` in `GAOptimize.m` with the desired function.
2. The selection function can be changed for both survivors and parents:
   - `selection_elitism(fitness, P);`
   - `selection_roulette(fitness, P);`
   - `selection_etournament(fitness, P, m);`
   - `selection_rtournament(fitness, P, m);`
   
   For survivors, replace the function on line `71` in `GAOptimize.m` with `selection_[NAME](tmp_fit, P[, m]);`.
   For parents, replace the functions on lines `85` and `87` in `GAOptimize.m` with `selection_[NAME](fitness, P[, m]);`.
   *Note:* `m` indicates the tournament size and must be replaced with a positive integer `≤ (P - S)/2`.
3. The crossover algorithm used to generate child genomes can be changes:
   - `crossover_1point(gen_p1, gen_p2, n);`
   - `crossover_2point(gen_p1, gen_p2, n);`
   - `crossover_probselect(gen_p1, gen_p2, n, m);`
   
   Replace the function on line `97` in `GAOptimize.m` with `crossover_[NAME](gen_p1, genP2, n[, m]);`.
   *Note:* `m` is `P(child(i) = gen_p1(i))`. Smaller values of `m` create a bias towards the first parent. Ensure `0 < m < 1`.
