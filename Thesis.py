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

def meiosis(arg,v,f):
    gametes = []
    if len(arg) == 2:  
        chromosomes = random.choices([1,2],weights=(1-v,v),k=1) # choice for reduced or unreduced gametes
        if chromosomes[0] == 1:
            gametes.append(random.sample(arg,1)[0])
        else:
            gametes.append(arg[0])
            gametes.append(arg[1])
    elif len(arg) == 3 and random.uniform(0,1) < f:
        chromosomes = random.choices([1,2,3],weights=(0.25,0.25,0.50),k=1) # 1 represents haploid, 2 diploid and 3 triploid (optional)
        poss_gametes = random.sample(arg,chromosomes[0])
        if chromosomes[0] == 1:
            gametes.append(poss_gametes[0])
        elif chromosomes[0] == 2:
            gametes.append(poss_gametes[0])
            gametes.append(poss_gametes[1])
        elif chromosomes[0] == 3:
            gametes.append(poss_gametes[0])
            gametes.append(poss_gametes[1])
            gametes.append(poss_gametes[2])
    elif len(arg) == 4:
        chromosomes = random.choices([0,2],weights=(v,1-v),k=1) # choice for diploid gametes or no gametes
        if chromosomes[0] == 2:
            poss_gametes = random.sample(arg,2)
            gametes.append(poss_gametes[0])
            gametes.append(poss_gametes[1])   
    return tuple(gametes)

def get_mean_of_freq_overreps(data,reps):
    mean_values = [0] * len(data[0])
    for rep in data:
        for i, value in enumerate(rep):
            mean_values[i] += value
    mean_values = [(value/reps) for value in mean_values]
    return mean_values

def recombination(individual, recombination_frequency = 0.3, length = 2):
    individual = sorted(individual, key=len, reverse=True) # sort individual so first chromosome set is always longest
    recombinated_individual = []
    # case of 1 x 1
    if (random.uniform(0,1) > recombination_frequency):
        random_index = random.randint(1,length -1)
        if len(individual[0]) == 1 and len(individual[1]) == 1:
            chromo_1 = individual[0]
            chromo_2 = individual[1]
            recombinated_chromo_1 = chromo_1[0][:random_index] + chromo_2[0][random_index:]
            recombinated_chromo_2 = chromo_2[0][:random_index] + chromo_1[0][random_index:]
            recombinated_individual.append(recombinated_chromo_1)
            recombinated_individual.append(recombinated_chromo_2)
        # case of 2 x 1
        if len(individual[0]) == 2 and len(individual[1]) == 1:
            chromo_1 = random.choice(individual[0]) # choose 1 chromosome from first set                
            remaining_chromo1 = [chromo for chromo in individual[0] if chromo != chromo_1][0]           
            chromo_2 = random.choice(individual[1]) # choose 1 chromsome from other set                 
            recombinated_chromo_1 = chromo_1[:random_index] + chromo_2[random_index:]               
            recombinated_chromo_2 = chromo_2[:random_index] + chromo_1[random_index:]               
            recombinated_individual.append([recombinated_chromo_1] + [remaining_chromo1])           
            recombinated_individual.append([recombinated_chromo_2])                                   
        # case of 2 x 2
        if len(individual[0]) == 2 and len(individual[1]) == 2:
            chromo_1 = random.choice(individual[0]) # choose 1 chromosome from first set                # ['A', 'A']
            remaining_chromo1 = [chromo for chromo in individual[0] if chromo != chromo_1][0]           # ['a', 'a']
            chromo_2 = random.choice(individual[1]) # choose 1 chromsome from other set                 # ['b', 'b']
            remaining_chromo2 = [chromo for chromo in individual[1] if chromo != chromo_2][0]           # ['B', 'B']
            recombinated_chromo_1 = chromo_1[:random_index] + chromo_2[random_index:]               # ['A', 'b']
            recombinated_chromo_2 = chromo_2[:random_index] + chromo_1[random_index:]               # ['b', 'A']
            recombinated_individual.append([recombinated_chromo_1] + [remaining_chromo1])           # [['A', 'b'],['a', 'a']]
            recombinated_individual.append([recombinated_chromo_2] + [remaining_chromo2])           # [['b', 'A'],['B', 'B']]
        # case of 3 x 1
        if len(individual[0]) == 3 and len(individual[1]) == 1:
            chromo_1 = random.choice(individual[0]) # choose 1 chromosome from first set                
            remaining_chromo1 = [chromo for chromo in individual[0] if chromo != chromo_1]           
            chromo_2 = random.choice(individual[1]) # choose 1 chromsome from other set                         
            recombinated_chromo_1 = chromo_1[:random_index] + chromo_2[random_index:]               
            recombinated_chromo_2 = chromo_2[:random_index] + chromo_1[random_index:]
            recombinated_individual.append([recombinated_chromo_1] + remaining_chromo1)
            recombinated_individual.append([recombinated_chromo_2])
    else: 
        recombinated_individual = individual
    return recombinated_individual

def cytotype_dynamics(initial_size:int,ploidy_number:int,loci_number:int,max_generations:int,v:float,f:float,reps:int):
    diploids = []
    triploids = []
    tetraploids = []
    for i in range(0,reps):
        initial_population = initialize_population(initial_size,loci_number,ploidy_number) # na each rep opnieuw initialpop opstarten ==> dus een nieuwe onafhankelijke trials
        t = 1
        n = len(initial_population)
        freq_diploid, freq_triploid, freq_tetraploid = [1], [0], [0]
        new_pop = initial_population
        while t < max_generations:
            offsprings = []
            while len(offsprings) < n:
                pair = random.sample(new_pop,2)
                first = meiosis(pair[0],v,f)
                second = meiosis(pair[1],v,f)
                if first != () and second != () and len(list(first + second)) < 5:
                    #offspring = recombination([first] + [second])
                    offspring = list(first + second)
                    offsprings.append(offspring)
            new_pop = offsprings
            freq_diploid.append(sum(1 for individual in new_pop if len(individual) == 2)/len(new_pop)) # count how many diploids in each generation
            freq_triploid.append(sum(1 for individual in new_pop if len(individual) == 3)/len(new_pop)) # count how many triploids in each generation
            freq_tetraploid.append(sum(1 for individual in new_pop if len(individual) == 4)/len(new_pop)) # count how many tetraploids in each generation
            t += 1
        diploids.append(freq_diploid)
        triploids.append(freq_triploid)
        tetraploids.append(freq_tetraploid)
    mean_diploids = get_mean_of_freq_overreps(diploids,reps)
    mean_triploids = get_mean_of_freq_overreps(triploids,reps)
    mean_tetraploids = get_mean_of_freq_overreps(tetraploids,reps)
    return max_generations,mean_diploids,mean_triploids,mean_tetraploids,diploids,triploids,tetraploids

def makeplot(data,data2):
    t = list(range(1,data[0]+1))
    plt.figure(figsize=(12, 5))

    plt.subplot(1,2,1)
    plt.plot(t,data[1], label= "Frequency diploid")
    plt.plot(t,data[2], label= "Frequency triploid")
    plt.plot(t,data[3], label= "Frequency tetraploid")
    plt.xlabel("Generation")
    plt.ylabel("Frequency")
    plt.title("Frequencies of di-,tri- and tetraploids per generation")
    plt.legend()

    plt.subplot(1,2,2)
    plt.plot(t,data2[1], label= "Frequency diploid")
    plt.plot(t,data2[2], label= "Frequency triploid")
    plt.plot(t,data2[3], label= "Frequency tetraploid")
    plt.xlabel("Generation")
    plt.ylabel("Frequency")
    plt.title("Frequencies of di-,tri- and tetraploids per generation")
    plt.legend()
    plt.show()

def recursive_NL_equations(max_generations,v,f,d1,d2,d3):
    # Initial conditions
    x2 = 1  # represent frequencies of diploid
    x3 = 0  # represent frequencies of triploid
    x4 = 0  # represent frequencies of tetraploid
    x2_list = [1]
    x3_list = [0]
    x4_list = [0]
    for i in range(max_generations-1):
        # update gamete proportions
        g1 =x2 * (1-v) + d1 * x3 * f
        g2= x2 * v + d2 * x3 * f + x4 * (1-v)
        g3 = d3 * x3 * f
        o = (g1 + g2 + g3)**2
        # update genotype frequencies
        x2 = g1**2 /o
        x3 = 2 * g1 * g2 / o
        x4 = (g2**2 + (2 * g1 * g3)) / o
        x2_list.append(x2)
        x3_list.append(x3)
        x4_list.append(x4)
    return max_generations,x2_list,x3_list,x4_list

makeplot(cytotype_dynamics(100,2,2,max_generations=100,v=0.01,f=0.3,reps=1000),recursive_NL_equations(max_generations=100,v=0.01,f=0.3,d1=0.25,d2=0.25,d3=0.5))


