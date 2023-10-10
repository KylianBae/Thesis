# Masterthesis
import random
import matplotlib.pyplot as plt

def initialize_population(initial_size:int,ploidy_number:int,loci_number:int):
    """_summary_
    Initialize the starting population with a given number of individuals of the population (initial_size), with its ploidy level(ploidy_number)
    and number of loci per chromosome (loci_number)
    Args:
        initial_size (int): Size of the initial population of individuals to start with
        ploidy_number (int): Number of chromosome sets per individual
        loci_number (int): Number of loci per chromosome

    Returns:
        List of lists which represents the population with its individuals based on the given parameters
    """
    population = []
    for i in range(0, initial_size):
        individual = []
        for i in range(0, ploidy_number):
            chromosome_set = []
            for i in range(0, loci_number):
                chromosome_set.append(random.choice([0,1]))
            individual.append(chromosome_set)
        population.append(individual)
    return population

def meiosis(arg,v):
    gametes = []
    if len(arg) == 2:  
        chromosomes = random.choices([1,2],weights=(1-v,v),k=1) # choice for reduced or unreduced gametes
        if chromosomes[0] == 1:
            gametes.append(random.sample(arg,1)[0])
        else:
            gametes.append(arg[0])
            gametes.append(arg[1])
    if len(arg) == 3 and random.uniform(0,1) < 0.3:
        chromosomes = random.choices([1,2,3],weights=(25,25,50),k=1) # 1 represents haploid, 2 diploid and 3 triploid (optional)
        poss_gametes = random.sample(arg,chromosomes[0])
        if chromosomes[0] == 1:
            gametes.append(poss_gametes[0])
        if chromosomes[0] == 2:
            gametes.append(poss_gametes[0])
            gametes.append(poss_gametes[1])
        if chromosomes[0] == 3:
            gametes.append(poss_gametes[0])
            gametes.append(poss_gametes[1])
            gametes.append(poss_gametes[2])
    if len(arg) == 4:
        chromosomes = random.choices([1,2],weights=(1- v,v),k=1) # choice for diploid gametes or no gametes
        if chromosomes[0] == 1:
            poss_gametes = random.sample(arg,2)
            gametes.append(poss_gametes[0])
            gametes.append(poss_gametes[1])   
    return gametes

def cytotype_dynamics(initial_population:list,max_generations:int,v:float, reps:int):
    results = []
    diploids = []
    triploids = []
    tetraploids = []
    for i in range(0,reps):
        t = 1
        n = len(initial_population)
        freq_diploid, freq_triploid, freq_tetraploid = [n], [0], [0]
        while t < max_generations:
            new_pop = []
            while len(new_pop) < n:
                pair = random.sample(initial_population,2)
                first = meiosis(pair[0],v)
                second = meiosis(pair[1],v)
                if first != [] and second != [] and len(first + second) < 5:
                    offspring = first + second
                    new_pop.append(offspring)
            freq_diploid.append(sum(1 for individual in new_pop if len(individual) == 2)) # count how many diploids in each generation
            freq_triploid.append(sum(1 for individual in new_pop if len(individual) == 3)) # count how many triploids in each generation
            freq_tetraploid.append(sum(1 for individual in new_pop if len(individual) == 4)) # count how many tetraploids in each generation
            t += 1
        diploids.append(freq_diploid)
        triploids.append(freq_triploid)
        tetraploids.append(freq_tetraploid)
    mean_diploids = get_mean_of_freq_overreps(diploids,reps)
    mean_triploids = get_mean_of_freq_overreps(triploids,reps)
    mean_tetraploids = get_mean_of_freq_overreps(tetraploids,reps)
    return max_generations,mean_diploids,mean_triploids,mean_tetraploids

def get_mean_of_freq_overreps(data,reps):
    mean_values = [0] * len(data[0])
    for rep in data:
        for i, value in enumerate(rep):
            mean_values[i] += value
    mean_values = [(value/reps) for value in mean_values]
    return mean_values


def makeplot(data):
    t = list(range(1,data[0]+1))
    plt.figure()
    plt.plot(t,data[1], label= "Frequency diploid")
    plt.plot(t,data[2], label= "Frequency triploid")
    plt.plot(t,data[3], label= "Frequency tetraploid")
    plt.xlabel("Generation")
    plt.ylabel("Frequency")
    plt.title("Plot of frequencies of diploids, triploids and tetraploids per generation")
    plt.legend()
    plt.show()

# Example
x = initialize_population(100,2,2)
print(makeplot(cytotype_dynamics(x,100,0.1,20)))
