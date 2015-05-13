A Simple Meta Genetic Algorithm
==========================
This a Meta-Genetic Algorithm about optimizing crossover and mutation probabilities in SGA presented [here](https://github.com/panvagenas/simple-genetic-algorithm)

Options
-------

**Usage:** meta_sga [options]

**Optional arguments:**

  `-h, --help`                                          show this help message and exit
  
  `-gp POPSIZE, --popsize POPSIZE`                      GA Population Size
  
  `-gg MAXGENS, --maxgens MAXGENS`                      GA Generations
  
  `-gn NVARS, --nvars NVARS`                            GA N Vars
  
  `-ga ANTIGENPOPSIZE, --antigenpopsize ANTIGENPOPSIZE` Antigens Population Size
  
  `-gr GA_RUNTIMES, --ga_runtimes GA_RUNTIMES`          Number of times GA will be called in order to eval
                                                        Meta-GA chromosome
                                                        
  `-g, --run_ga`                                        Run Only SGA
  
  `-mp META_POPSIZE, --meta_popsize META_POPSIZE`       Meta-GA Population Size
  
  `-mg META_MAXGENS, --meta_maxgens META_MAXGENS`       Meta-GA Generations
  
  `-mn META_NVARS, --meta_nvars META_NVARS`             Meta-GA N Vars
  
  `-px CROSSOVER, --crossover CROSSOVER`                Meta-GA Crossover Probability
  
  `-pm MUTATION, --mutation MUTATION`                   Meta-GA Mutation Probability
  
  `-v, --verbose`                                       Verbose info to std output
  
  `-sp SPLIT_POINT, --split_point SPLIT_POINT`          Meta-GA Chromosome split point
  
  `-e, --elitism`                                       Use elitism (effects only Meta-GA)
 
*this is a work in progress*