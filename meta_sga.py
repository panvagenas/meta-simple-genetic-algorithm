#!/usr/bin/python3
import argparse
from abc import ABCMeta, abstractmethod
from random import random, randrange, uniform, randint
from math import ceil
from copy import copy, deepcopy
from time import time


parser = argparse.ArgumentParser(description='Meta-SGA Algorithm')

parser.add_argument('-gp', '--popsize', dest='popsize', default=50, type=int,
                    help='GA Population Size')
parser.add_argument('-gg', '--maxgens', dest='maxgens', default=20, type=int,
                    help='GA Generations')
parser.add_argument('-gn', '--nvars', dest='nvars', default=8, type=int,
                    help='GA N Vars')
parser.add_argument('-ga', '--antigenpopsize', dest='antigenpopsize', default=1000, type=int,
                    help='Antigens Population Size')
parser.add_argument('-gr', '--ga_runtimes', dest='ga_runtimes', default=1, type=int,
                    help='Number of times GA will be called in order to eval Meta-GA chromosome')
parser.add_argument('-g', '--run_ga', dest='run_ga', action='store_const', const=True,
                    help='Run Only SGA')

parser.add_argument('-mp', '--meta_popsize', dest='meta_popsize', default=100, type=int,
                    help='Meta-GA Population Size')
parser.add_argument('-mg', '--meta_maxgens', dest='meta_maxgens', default=50, type=int,
                    help='Meta-GA Generations')
parser.add_argument('-mn', '--meta_nvars', dest='meta_nvars', default=8, type=int,
                    help='Meta-GA N Vars')

parser.add_argument('-px', '--crossover', dest='crossover', default=0.7, type=float,
                    help='Meta-GA Crossover Probability')
parser.add_argument('-pm', '--mutation', dest='mutation', default=0.01, type=float,
                    help='Meta-GA Mutation Probability')
parser.add_argument('-v', '--verbose', dest='verbose', action='store_const', const=True,
                    help='Verbose info to std output')
parser.add_argument('-sp', '--split_point', dest='split_point', default=0, type=int,
                    help='Meta-GA Chromosome split point')
parser.add_argument('-e', '--elitism', dest='elitism', action='store_const', const=True,
                    help='Use elitism (effects only Meta-GA)')

args = parser.parse_args()


class Chromosome:
    def __init__(self):
        self.fitness = 0
        self.gene = ''
        self.mate = False
        self.selection_prob = 0
        self.cumulative_prob = 0


# Abstract population class
# Meta Population and GA Population use this as base class
class Population(object):
    __metaclass__ = ABCMeta

    def __init__(self, nvars, pop_size):
        self.best = Chromosome()
        self.population = [Chromosome() for x in range(pop_size)]
        self.nvars = nvars
        self.pop_size = pop_size
        self.history = []
        self.mean_history = []

        for chromosome in self.population:
            chromosome.gene = "".join([str(randint(0, 1)) for j in range(0, nvars)])

    @abstractmethod
    def evaluate(self):
        pass

    @abstractmethod
    def crossover(self, pxover):
        pass

    def select_new(self):
        buffered_pop = self.population[:]

        total_fitness = 0
        for member in self.population:
            total_fitness += member.fitness

        for member in self.population:
            member.selection_prob = member.fitness / total_fitness

        self.population[0].cumulative_prob = self.population[0].selection_prob

        for i in range(1, self.pop_size):
            self.population[i].cumulative_prob = self.population[i - 1].cumulative_prob + self.population[
                i].selection_prob

        for spin_num in range(0, self.pop_size):
            roulette = random()

            if roulette <= self.population[0].cumulative_prob:
                buffered_pop[spin_num] = copy(self.population[0])
            else:
                for mem in range(1, self.pop_size):
                    if self.population[i - 1].cumulative_prob < roulette <= self.population[i].cumulative_prob:
                        buffered_pop[spin_num] = copy(self.population[i])
                        break

        self.population = buffered_pop

    def mutate(self, pmutation):
        for member in self.population:
            gene = list(member.gene)
            for i in range(0, len(gene)):
                if random() < pmutation:
                    if gene[i] == '0':
                        gene[i] = '1'
                    else:
                        gene[i] = '0'
            member.gene = "".join(gene)


class GAPopulation(Population):
    def __init__(self, nvars, pop_size, antigene_pop_size):
        Population.__init__(self, nvars, pop_size)
        self.antigene_pop_size = antigene_pop_size
        self.antigene_dna = ''
        self.evaluation_dict = {}

        self.antigene_dna = "".join([str(randint(0, 1)) for j in range(0, nvars * self.antigene_pop_size)])

    def evaluate(self):
        mean = 0

        for mem in self.population:
            if mem.gene in self.evaluation_dict:
                mem.fitness = self.evaluation_dict.get(mem.gene)
            else:
                mem_dna = mem.gene * self.antigene_pop_size

                xor = int(mem_dna, 2) ^ int(self.antigene_dna, 2)
                mem.fitness = bin(xor).count("1")

                self.evaluation_dict[mem.gene] = mem.fitness

                if mem.fitness > self.best.fitness:
                    self.best = copy(mem)

            mean += mem.fitness

        self.history.append(self.best.fitness)
        self.mean_history.append(mean / self.pop_size)

        return self.best.fitness

    def crossover(self, pxover):
        lovers = 0

        for chromosome in self.population:
            if random() < pxover:
                chromosome.mate = True
                lovers += 1
            else:
                chromosome.mate = False

        if lovers == self.pop_size and lovers % 2 != 0:
            lovers -= 1

        if lovers % 2 != 0:
            while True:
                chromosome = self.population[randrange(0, self.pop_size)]
                if not chromosome.mate:
                    chromosome.mate = True
                    lovers += 1
                    break

        parent_1 = 0
        couples = int(ceil(lovers / 2))
        for c in range(0, couples):
            while not self.population[parent_1].mate:
                parent_1 += 1

            parent_2 = parent_1 + 1
            while not self.population[parent_2].mate:
                parent_2 += 1

            xover_point = randrange(1, self.nvars)

            gene_1 = self.population[parent_1].gene[:]
            gene_2 = self.population[parent_2].gene[:]

            self.population[parent_1].gene = "".join([gene_1[:xover_point], gene_2[xover_point:]])
            self.population[parent_2].gene = "".join([gene_2[:xover_point], gene_1[xover_point:]])

            parent_1 = parent_2 + 1


class MetaPopulation(Population):
    def __init__(self, nvars, pop_size, ga_nvars, ga_pop_size, ga_antigene_pop_size, ga_maxgens, split_point):
        Population.__init__(self, nvars, pop_size)
        self.ga_maxgens = ga_maxgens
        self.ga_nvars = ga_nvars
        self.ga_pop_size = ga_pop_size
        self.ga_antigene_pop_size = ga_antigene_pop_size
        self.split_point = split_point
        self.px_history = []
        self.pm_history = []

    def evaluate(self):
        mean = 0

        for mem in self.population:
            mem.fitness = 0

            if not args.elitism:
                self.best = Chromosome()

            mean_for_px_pm = 0.0
            px = self.px_from_chromosome(mem)
            pm = self.pm_from_chromosome(mem)

            for I in range(0, args.ga_runtimes):
                ga_pop = self.run_ga(
                    px, pm,
                    self.ga_pop_size,
                    self.ga_antigene_pop_size,
                    self.ga_nvars,
                    self.ga_maxgens
                )
                mean_for_px_pm += ga_pop.best.fitness

            mem.fitness = mean_for_px_pm / args.ga_runtimes

            if mem.fitness > self.best.fitness:
                self.best = copy(mem)

            mean += mem.fitness

        self.history.append(self.best.fitness)
        self.mean_history.append(mean / self.pop_size)
        self.px_history.append(self.px_from_chromosome(self.best))
        self.pm_history.append(self.pm_from_chromosome(self.best))

        return self.best.fitness

    def crossover(self, pxover):
        lovers = 0

        for chromosome in self.population:
            if random() < pxover:
                chromosome.mate = True
                lovers += 1
            else:
                chromosome.mate = False

        if lovers == self.pop_size and lovers % 2 != 0:
            lovers -= 1

        if lovers % 2 != 0:
            while True:
                chromosome = self.population[randrange(0, self.pop_size)]
                if not chromosome.mate:
                    chromosome.mate = True
                    lovers += 1
                    break

        parent_1 = 0
        couples = int(ceil(lovers / 2))
        for c in range(0, couples):
            while not self.population[parent_1].mate:
                parent_1 += 1

            parent_2 = parent_1 + 1
            while not self.population[parent_2].mate:
                parent_2 += 1

            px_xover_point = randrange(1, self.split_point)
            pm_xover_point = randrange(1, self.nvars - self.split_point)

            px_gene_1 = self.population[parent_1].gene[:self.split_point]
            px_gene_2 = self.population[parent_2].gene[:self.split_point]
            pm_gene_1 = self.population[parent_1].gene[self.split_point:]
            pm_gene_2 = self.population[parent_2].gene[self.split_point:]

            new_gene_1 = "".join([
                px_gene_1[:px_xover_point],
                px_gene_2[px_xover_point:],
                pm_gene_1[:pm_xover_point],
                pm_gene_2[pm_xover_point:]
            ])

            new_gene_2 = "".join([
                px_gene_2[:px_xover_point],
                px_gene_1[px_xover_point:],
                pm_gene_2[:pm_xover_point],
                pm_gene_1[pm_xover_point:]
            ])

            self.population[parent_1].gene = new_gene_1
            self.population[parent_2].gene = new_gene_2

            parent_1 = parent_2 + 1

    def px_from_chromosome(self, chromosome):
        px_bin = chromosome.gene[:self.split_point]
        px = 0.0
        exp = 1
        for value in list(px_bin):
            px += (1 / (2 ** exp)) * int(value)
            exp += 1

        return px

    def pm_from_chromosome(self, chromosome):
        pm_bin = chromosome.gene[self.split_point:]
        pm = 0.0
        exp = 1
        for value in list(pm_bin):
            pm += (1 / (2 ** exp)) * int(value)
            exp += 1
        return pm

    @staticmethod
    def run_ga(px, pm, pop_size, antigene_pop_size, nvars, maxgens):
        pop = GAPopulation(nvars, pop_size, antigene_pop_size)

        ev = pop.evaluate()

        for generation in range(0, maxgens):
            pop.select_new()
            pop.crossover(px)
            pop.mutate(pm)
            pop.evaluate()

        return pop


def run_sga():
    print('Running SGA Algorithm')

    mean = 0

    for i in range(0, args.ga_runtimes):
        best_pop = GAPopulation(args.nvars, args.popsize, args.antigenpopsize)

        pop = GAPopulation(args.nvars, args.popsize, args.antigenpopsize)
        pop.evaluate()

        for generation in range(0, args.maxgens):
            pop.select_new()
            pop.crossover(args.crossover)
            pop.mutate(args.mutation)
            pop.evaluate()

        mean += pop.best.fitness
        if best_pop.best.fitness < pop.best.fitness:
            best_pop = deepcopy(pop)

    print('Mean fitness: %f' % (mean / args.ga_runtimes))
    pass


def run_meta_sga():
    start = time()
    print('Running Meta-SGA')

    if not (0 < args.split_point < args.meta_nvars):
        args.split_point = int(ceil(args.meta_nvars / 2))

    pop = MetaPopulation(args.meta_nvars, args.meta_popsize, args.nvars, args.popsize, args.antigenpopsize,
                         args.popsize, args.split_point)
    pop.evaluate()

    for generation in range(0, args.meta_maxgens):
        print('Generation %d:' % generation)
        print('Best gene: %s' % pop.best.gene)
        print('Fitness: %f' % pop.best.fitness)
        print('Crossover probability: %f' % pop.px_from_chromosome(pop.best))
        print('Mutation probability: %f' % pop.pm_from_chromosome(pop.best))
        if args.verbose:
            print_px_pm_spread(pop)
        elapsed_time = time() - start
        print('Time taken: %f' % elapsed_time)
        time_remaining = elapsed_time * (args.meta_maxgens - generation)
        print('Estimated time remaining: %.1f min' % (time_remaining / 60))

        start = time()
        print('')
        print('Processing Meta Generation %d of %d' % (generation + 1, args.meta_maxgens))
        pop.select_new()
        pop.crossover(args.crossover)
        pop.mutate(args.mutation)
        pop.evaluate()

    print('Best Gene: %s' % pop.best.gene)
    print('For crossover probability: %f and mutation probability: %f' %
          (pop.px_from_chromosome(pop.best), pop.pm_from_chromosome(pop.best)))
    print('Fitness: %f' % pop.best.fitness)

    return pop


def print_px_pm_spread(population):
    spread_px = {i: 0 for i in range(0, 10)}
    spread_pm = {i: 0 for i in range(0, 10)}
    for mem in population.population:
        if population.px_from_chromosome(mem) > 0.9:
            spread_px[9] += 1
        elif population.px_from_chromosome(mem) > 0.8:
            spread_px[8] += 1
        elif population.px_from_chromosome(mem) > 0.7:
            spread_px[7] += 1
        elif population.px_from_chromosome(mem) > 0.6:
            spread_px[6] += 1
        elif population.px_from_chromosome(mem) > 0.5:
            spread_px[5] += 1
        elif population.px_from_chromosome(mem) > 0.4:
            spread_px[4] += 1
        elif population.px_from_chromosome(mem) > 0.3:
            spread_px[3] += 1
        elif population.px_from_chromosome(mem) > 0.2:
            spread_px[2] += 1
        elif population.px_from_chromosome(mem) > 0.1:
            spread_px[1] += 1
        else:
            spread_px[0] += 1

        if population.pm_from_chromosome(mem) > 0.9:
            spread_pm[9] += 1
        elif population.pm_from_chromosome(mem) > 0.8:
            spread_pm[8] += 1
        elif population.pm_from_chromosome(mem) > 0.7:
            spread_pm[7] += 1
        elif population.pm_from_chromosome(mem) > 0.6:
            spread_pm[6] += 1
        elif population.pm_from_chromosome(mem) > 0.5:
            spread_pm[5] += 1
        elif population.pm_from_chromosome(mem) > 0.4:
            spread_pm[4] += 1
        elif population.pm_from_chromosome(mem) > 0.3:
            spread_pm[3] += 1
        elif population.pm_from_chromosome(mem) > 0.2:
            spread_pm[2] += 1
        elif population.pm_from_chromosome(mem) > 0.1:
            spread_pm[1] += 1
        else:
            spread_pm[0] += 1

    print('Crossover Probability Spread: ')
    for i in spread_px:
        print('%d-%d\t%s' % (i, i + 1, (spread_px[i] * '=')))
    print('Mutation Probability Spread: ')
    for i in spread_px:
        print('%d-%d\t%s' % (i, i + 1, (spread_pm[i] * '=')))


if args.run_ga:
    run_sga()
else:
    pop = run_meta_sga()

print('\n')
print('=' * 19, ' Arguments ', '=' * 19)
print(
    ('GA Population Size: %d\n'
     'GA Generations: %d\n'
     'GA N Vars: %d\n'
     'GA Antigens Population Size: %d\n'
     'GA Runtimes: %d\n'
     'Crossover Probability: %f\n'
     'Mutation Probability: %f\n' %
     (args.popsize, args.maxgens, args.nvars, args.antigenpopsize, args.ga_runtimes, args.crossover, args.mutation))
)

if not args.run_ga:
    if args.elitism:
        elitism = 'Yes'
    else:
        elitism = 'No'
    print('Meta-GA Population Size: %d\n'
          'Meta-GA Generations: %d\n'
          'Meta-GA N Vars: %d\n'
          'Meta-GA Chromosome split point: %d\n'
          'Using Elitism: %s\n' %
          (args.meta_popsize, args.meta_maxgens, args.meta_nvars, args.split_point, elitism))

    print('=' * 20, ' Results ', '=' * 20)
    print('Best Gene: %s' % pop.best.gene)
    print('For crossover probability: %f and mutation probability: %f' %
          (pop.px_from_chromosome(pop.best), pop.pm_from_chromosome(pop.best)))
    print('Fitness: %f' % pop.best.fitness)

print('=' * 51, '\n')