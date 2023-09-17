# Masterthesis
import random
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

def meiosis(individuals,v=0.5):
    gametes = []
    for arg in individuals:
        if len(arg) == 2:
            chromosomes = random.choices([1,2],weights=(1-v,v),k=1)
            if chromosomes[0] == 1:
                gametes.append(random.sample(arg,1)[0])
            else:
                gametes.append(arg)
        if len(arg) == 3 and random.uniform(0,1) < 0.3:
            chromosomes = random.choices([1,2,3],weights=(25,25,50),k=1) # 1 represents haploid, 2 diploid and 3 triploid (optional)
            gametes.append(random.sample(arg,chromosomes[0])[0])
        if len(arg) == 4:
            chromosomes = random.choices([1,2],weights=(1- v,v),k=1)
            if chromosomes[0] == 1:
                gametes.append(random.sample(arg,1)[0])    
    return flattenlist(gametes)

def flattenlist(lst):
    flat = []
    for item in lst:
        if isinstance(item,list):
            flat.extend(flattenlist(item))
        else:
            flat.append(item)
    return flat

def mate(start_population):
    new_pop = []
    while len(new_pop) < len(start_population):
        pair = random.sample(start_population,2)
        gametes = meiosis(pair)
        if len(gametes) >= 2:
           new_pop.append(gametes)
    return new_pop

x = initialize_population(4,2,2)
print(x)
print(mate(x))
