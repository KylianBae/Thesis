# Masterthesis
import random
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

def recursive_NL_equations(max_generations,v,f,d1,d2,d3):
    """ Within this function the mathematical model is run using equations 1-7 (see manuscript) which will simulate the cytotype dynamics
    within an infinite size population for a determined amount of generations using the given parameters v,f,d1,d2,d3

    Args:
        max_generations (int): maximum amount of generations the model has to run for
        v (float): probability of unreduced gamete production 
        f (float): triploid fitness
        d1 (float): haploid gamete production from triploid individual
        d2 (float): diploid gamete production from triploid individual
        d3 (float): triploid gamete production from triploid individual

    Returns:
        list: list of cytotype frequencies for each generation within the mathematical model
    """
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

def plot_line(inp,title):
    """ Function to plot the cytotype frequencies per generation

    Args:
        inp (list): cytotype frequencies for each generation
        title (str): title of plot
    """
    x_values =  list(range(1,inp[0]+1))
    diploid_values = inp[1]
    triploid_values = inp[2]
    tetraploid_values = inp[3]
    # Create a line plot using seaborn
    sns.set(style="whitegrid")  # Set the style of the plot

    # Create a line plot
    sns.lineplot(x=x_values, y=diploid_values)
    sns.lineplot(x=x_values, y=triploid_values)
    sns.lineplot(x=x_values, y=tetraploid_values)

    # Set labels and title
    plt.xlabel("Generation")
    plt.ylabel("Frequency")
    plt.title(title)
    plt.legend()

    # Show the plot
    plt.show()

def plot_fertility_unreduced_gametes():
    """Function to plot the cytotype dynamics with the unreduced gamete production on x-axis and equilibrium frequency at y-axis per triploid fitness
    """
    def freq_fertility(fertility):
        max_generations = 100
        v_rates = [i * 0.001 for i in range(0, 251)]
        Eq_freq_2,Eq_freq_3,Eq_freq_4 = [],[],[]
        for i in v_rates:
            #x = recursive_NL_equations(max_generations=max_generations,v=i,f=fertility,d1=0.25,d2=0.25,d3=0.5) # with triploid gametes
            x = recursive_NL_equations(max_generations=max_generations,v=i,f=fertility,d1=0.5,d2=0.5,d3=0)  # without triploid gametes
            Eq_freq_2.append(x[1][-1])
            Eq_freq_3.append(x[2][-1])
            Eq_freq_4.append(x[3][-1])
        return [Eq_freq_2,Eq_freq_3,Eq_freq_4]

    def plot_line(x_values,inp,inp1,inp2,inp3,inp4):
        diploid_values = inp[0]
        triploid_values = inp[1]
        tetraploid_values = inp[2]

        diploid_values1 = inp1[0]
        triploid_values1 = inp1[1]
        tetraploid_values1 = inp1[2]

        diploid_values2 = inp2[0]
        triploid_values2 = inp2[1]
        tetraploid_values2 = inp2[2]

        diploid_values3 = inp3[0]
        triploid_values3 = inp3[1]
        tetraploid_values3 = inp3[2]

        diploid_values4 = inp4[0]
        triploid_values4 = inp4[1]
        tetraploid_values4 = inp4[2]
        # Create a line plot using seaborn
        sns.set(style="whitegrid")  # Set the style of the plot

        # Create a line plot
        sns.lineplot(x=x_values, y=diploid_values, color='red')
        sns.lineplot(x=x_values, y=triploid_values, color='green')
        sns.lineplot(x=x_values, y=tetraploid_values, color='blue')

        sns.lineplot(x=x_values, y=diploid_values1, color='red')
        sns.lineplot(x=x_values, y=triploid_values1, color='green')
        sns.lineplot(x=x_values, y=tetraploid_values1,color='blue')

        sns.lineplot(x=x_values, y=diploid_values2, color='red')
        sns.lineplot(x=x_values, y=triploid_values2, color='green')
        sns.lineplot(x=x_values, y=tetraploid_values2,color='blue')

        sns.lineplot(x=x_values, y=diploid_values3, color='red')
        sns.lineplot(x=x_values, y=triploid_values3, color='green')
        sns.lineplot(x=x_values, y=tetraploid_values3,color='blue')

        sns.lineplot(x=x_values, y=diploid_values4, color='red')
        sns.lineplot(x=x_values, y=triploid_values4, color='green')
        sns.lineplot(x=x_values, y=tetraploid_values4,color='blue')

        # Set labels and title
        plt.xlabel("Probability of unreduced gametes v")
        plt.ylabel(" Equilibrium Frequency")
        plt.yticks([0, 0.25, 0.5, 0.75, 1])
        plt.legend()

        # Show the plot
        plt.show()

    v_rates = [i * 0.001 for i in range(0, 251)]
    inp = freq_fertility(0)
    inp1 = freq_fertility(0.1)
    inp2 = freq_fertility(0.2)
    inp3 = freq_fertility(0.3)
    inp4 = freq_fertility(0.4)

    plot_line(v_rates,inp,inp1,inp2,inp3,inp4)

def cytotypespecific_frequencies_di_tri_tetraploid(data_simulation,data_equations):
    """Plots the cytotype frequencies per generation per repetition for each cytotype separately
    Args:
        data_simulation (lst): list containing the cytotype frequencies simulated by the IBM
        data_equations (lst): list containing the cytotype frequencies simulated by the deterministic model
    """
    t = list(range(1,data_simulation[0]+1))
    t2 = list(range(1,data_equations[0]+1))

    fig,axes = plt.subplots(3,1,sharex=True)
    fig.suptitle("Cytotype-specific dynamics")

    # Diploid
    for rep in data_simulation[4]:
        sns.scatterplot(ax= axes[0],x=t, y=rep, color="red")
    sns.lineplot(ax= axes[0],x=t2, y=data_equations[1], color= "black")
    plt.ylim(0.6, 1)

    # Triploid
    for rep in data_simulation[5]:
        sns.scatterplot(ax= axes[1], x=t, y=rep, color="green")
    sns.lineplot(ax= axes[1], x=t2, y=data_equations[2], color= "black")
    plt.ylim(0, 0.4)

    # Tetraploid 
    for rep in data_simulation[6]:
        sns.scatterplot(ax= axes[2], x=t, y=rep, color="blue")
    sns.lineplot(ax= axes[2], x=t2, y=data_equations[3], color= "black")
    plt.ylim(0, 0.2)

    plt.show()

def phase_portrait(f= 0.3,v1=0.01,v3=0.07,v5=0.05,v10=0.14):
    """Function the form the phase portrait to study the dynamics of the haploid gamete frequency, bifurcation analysis

    Args:
        f (float, optional): triploid fitness. Defaults to 0.3.
        v1 (float, optional): first probability of unreduced gamete production. Defaults to 0.01.
        v3 (float, optional): second probability of unreduced gamete production. Defaults to 0.07.
        v5 (float, optional): third probability of unreduced gamete production. Defaults to 0.05.
        v10 (float, optional): fourth probability of unreduced gamete production. Defaults to 0.14.
    """
    # haploid gamete frequency: g^2(-2+v+2f) + g(+3-3v-3f) + (-1+v+f) = 0                
    # see manuscript 
    x = np.linspace(0,1.5,100)
    y1 = x**2 * (-2 + v1 + 2 * f) + x * (-3 * v1 - 3 * f + 3) + (-1 + v1 + f) 
    y3 = x**2 * (-2 + v3 + 2 * f) + x * (-3 * v3 - 3 * f + 3) + (-1 + v3 + f)
    #y5 = x**2 * (-2 + v5 + 2 * f) + x * (-3 * v5 - 3 * f + 3) + (-1 + v5 + f)
    y10 = x**2 * (-2 + v10 + 2 * f) + x * (-3 * v10 - 3 * f + 3) + (-1 + v10 + f)
    x_axis = 0 * x
    # Create a line plot using seaborn
    sns.set(style="whitegrid")  # Set the style of the plot

    # Create a line plot
    sns.lineplot(x=x, y=y1,label=f"v = {v1}")
    sns.lineplot(x=x, y=y3,label=f"v = {v3}")
    #sns.lineplot(x=x, y=y5,label=f"v = {v5}")
    sns.lineplot(x=x, y=y10,label=f"v = {v10}")
    sns.lineplot(x=x, y=x_axis, color="black")

    # Set labels and title
    plt.xlabel("Haploid gamete frequency")
    plt.ylabel("Rate of haploid gamete frequency change")
    plt.title("Phase portrait of haploid gamete frequency")
    plt.legend()
    #plt.ylim(-0.1, 0.1)
    # Show the plot
    plt.show()

def initialize_population(initial_size:int,ploidy_number:int,loci_number:int,type=None):
    """_summary_
    Initialize the starting population with a given number of individuals of the population (initial_size), with its ploidy level(ploidy_number)
    and number of loci per chromosome (loci_number). In this case 2 alleles exist for each loci. The population will be started to exist completely 
    of heterzygote individuals.
    Args:
        initial_size (int): Size of the initial population of individuals to start with
        ploidy_number (int): Number of chromosome sets per individual
        loci_number (int): Number of loci per chromosome
        type (int): type of initial population, None= population with different genotypes and fitnesses chosen at random, 0=population with worst possible fitness
        1= population with best possible fitness, 2=population with only heterozygous individuals, 3= population with only homozygous individuals, 4= population with only tetraploids of worst genotype
        5= population of triploids with worst genotype, 6= population mixed with genotypes of neutral fitness (O.8)


    Returns:
        List of lists which represents the population with its individuals based on the given parameters
    """
    population = []
    Heterozygous_possibilities = [[["A","B"],["a","b"]],[["a","b"],["A","B"]],[["A","b"],["a","B"]],[["a","B"],["A","b"]]]
    Homozygous_possibilities = [[["A","B"],["A","B"]],[["A","b"],["A","b"]],[["a","b"],["a","b"]],[["a","B"],["a","B"]]]
    if type == 0:
        while len(population) < initial_size:
            population.append([["a","B"],["a","B"]])
    if type == 1:
        while len(population) < initial_size:
            population.append([["A","b"],["A","b"]]) 
    if type == 2:
        while len(population) < initial_size:
            population.append(random.choice(Heterozygous_possibilities))
    if type == 3:
        while len(population) < initial_size:
            population.append(random.choice(Homozygous_possibilities))
    if type == 4:
        while len(population) < initial_size:
            population.append([["a","B"],["a","B"],["a","B"],["a","B"]])
    if type == 5:
        while len(population) < initial_size:
            population.append([["a","B"],["a","B"],["a","B"]])
    if type == 6:
        while len(population) < initial_size:
            population.append(random.choice([[["a","b"],["a","b"]],[["A","b"],["a","B"]],[["A","B"],["A","B"]]]))
    if type == None:
        while len(population) < initial_size:
            individual = []
            for i in range(0, ploidy_number):
                chromosome_set = []
                chromosome_set.append(random.choice(["A","a"]))
                chromosome_set.append(random.choice(["B","b"]))
                individual.append(chromosome_set)
            population.append(individual)
    return population

def meiosis(arg,v,f,recombination_frequency):
    """Simulated the process of meiosis in the IBM (see manuscript)

    Args:
        arg (lst): amount of individuals
        v (float): probability of unreduced gamete production
        recombination_frequency (float): rate at which recombination happens

    Returns:
        individual after meiosis
    """
    gametes = []
    arg = recombination(arg,recombination_frequency)
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
    """Function to get the mean cytotype frequencies per generation from all repetitions

    Args:
        data (list): list with equilibrium frequencies
        reps (int): amount of repetitions or runs of the IBM

    Returns:
        list containing mean cytotype frequencies per generation from all reps
    """
    mean_values = [0] * len(data[0])
    for rep in data:
        for i, value in enumerate(rep):
            mean_values[i] += value
    mean_values = [(value/reps) for value in mean_values]
    return mean_values

def get_mean_of_freq_overreps_counts(data,reps):
    """Function to get the mean cytotype frequencies per generation from all repetitions

    Args:
        data (list): list with equilibrium frequencies
        reps (int): amount of repetitions or runs of the IBM

    Returns:
        list containing mean cytotype frequencies per generation from all reps
    """
    mean_values = [0] * len(data[0])
    for rep in data:
        for i, value in enumerate(rep):
            if not isinstance(value,list):
                mean_values[i] += value
    mean_values = [(value/reps) for value in mean_values]
    return mean_values

def recombination(individual, recombination_frequency = 0.01, length = 2):
    """Implementation of the process of recombination (see manuscript)

    Args:
        individual (lst): _individual in which recombination should happen
        recombination_frequency (float, optional): probability that recombination happens per loci. Defaults to 0.01.
        length (int, optional): amount of loci. Defaults to 2.

    Returns:
        list: new individual with recombinated genetic information
    """
    individual = sorted(individual, key=len, reverse=True) # sort individual so first chromosome set is always longest
    # case of 2 gametes
    if len(individual) == 2:
        for i in range(0,length):
            if (random.uniform(0,1) < recombination_frequency):
                individual[0][i], individual[1][i] = individual[1][i], individual[0][i] 
    # case of 3 gametes
    if len(individual) == 3:
        for i in range(0,length):
            if (random.uniform(0,1) < recombination_frequency):
                individual[0][i], individual[1][i] = individual[1][i], individual[0][i]
        for i in range(0,length):
            if (random.uniform(0,1) < recombination_frequency):
                individual[1][i], individual[2][i] = individual[2][i], individual[1][i]
        for i in range(0,length):
            if (random.uniform(0,1) < recombination_frequency):
                individual[2][i], individual[0][i] = individual[0][i], individual[2][i]
    # case of 4 gametes
    if len(individual) == 4:
        for i in range(0,length):
            if (random.uniform(0,1) < recombination_frequency):
                individual[0][i], individual[1][i] = individual[1][i], individual[0][i]
        for i in range(0,length):
            if (random.uniform(0,1) < recombination_frequency):
                individual[1][i], individual[2][i] = individual[2][i], individual[1][i]
        for i in range(0,length):
            if (random.uniform(0,1) < recombination_frequency):
                individual[2][i], individual[0][i] = individual[0][i], individual[2][i]
        for i in range(0,length):
            if (random.uniform(0,1) < recombination_frequency):
                individual[3][i], individual[1][i] = individual[1][i], individual[3][i]
        for i in range(0,length):
            if (random.uniform(0,1) < recombination_frequency):
                individual[3][i], individual[2][i] = individual[2][i], individual[3][i]
        for i in range(0,length):
            if (random.uniform(0,1) < recombination_frequency):
                individual[3][i], individual[0][i] = individual[0][i], individual[3][i]
    recombinated_individual = individual
    return recombinated_individual

def cytotype_dynamics(initial_size:int,ploidy_number:int,loci_number:int,max_generations:int,v:float,f:float,reps:int,apply_fitness=True,selection_coeff=0.2,reference_fitness=0.8, recombination_frequency=0.01,type=None):
    """The individual based stochastic simulation as described in the manuscript

    Args:
        initial_size (int): initial population size (N)
        ploidy_number (int): ploidy level of initial population
        loci_number (int): amount of loci of each chromosome
        max_generations (int): maximum amount of generations the IBM should be run
        v (float): probability of unreduced gamete production
        f (float): triploid fitness
        reps (int): amount of repetitions the IBM should be run for
        apply_fitness (bool, optional): apply fitness function or not. Defaults to True.
        selection_coeff (float, optional): selection coefficient of fitness. Defaults to 0.2.
        reference_fitness (float, optional): reference fitness score when fitness is applied. Defaults to 0.8.
        recombination_frequency (float, optional): probability that recombination takes place per loci. Defaults to 0.01.
        type (_type_, optional): type of population. Defaults to None.

    Returns:
        mean cytotype frequencies and population per generation 
    """
    diploids = []
    triploids = []
    tetraploids = []
    population_per_generation_per_rep = []
    for i in range(0,reps):
        initial_population = initialize_population(initial_size,loci_number,ploidy_number,type=type) # na each rep opnieuw initialpop opstarten ==> dus een nieuwe onafhankelijke trials
        t = 1
        n = len(initial_population)
        freq_diploid, freq_triploid, freq_tetraploid = [1], [0], [0]
        new_pop = initial_population
        population_per_generation = [initial_population]
        while t < max_generations:
            offsprings = []
            while len(offsprings) < n:
                pair = random.sample(new_pop,2)
                if apply_fitness:
                    if fitness(pair[0],pair[1],selection_coeff=0.2,reference_fitness=0.8) == True:
                        first = meiosis(pair[0],v,f,recombination_frequency)
                        second = meiosis(pair[1],v,f,recombination_frequency)
                        if first != () and second != () and len(list(first + second)) < 5:
                            offspring = list(first + second)
                            offsprings.append(offspring)
                else:
                    first = meiosis(pair[0],v,f,recombination_frequency)
                    second = meiosis(pair[1],v,f,recombination_frequency)
                    if first != () and second != () and len(list(first + second)) < 5:
                        #offspring = recombination([first] + [second])
                        offspring = list(first + second)
                        offsprings.append(offspring)                  
            new_pop = offsprings
            population_per_generation.append(offsprings)
            freq_diploid.append(sum(1 for individual in new_pop if len(individual) == 2)/len(new_pop)) # count how many diploids in each generation
            freq_triploid.append(sum(1 for individual in new_pop if len(individual) == 3)/len(new_pop)) # count how many triploids in each generation
            freq_tetraploid.append(sum(1 for individual in new_pop if len(individual) == 4)/len(new_pop)) # count how many tetraploids in each generation
            t += 1
        diploids.append(freq_diploid)
        triploids.append(freq_triploid)
        tetraploids.append(freq_tetraploid)
        population_per_generation_per_rep.append(population_per_generation)
    mean_diploids = get_mean_of_freq_overreps(diploids,reps)
    mean_triploids = get_mean_of_freq_overreps(triploids,reps)
    mean_tetraploids = get_mean_of_freq_overreps(tetraploids,reps)
    return max_generations,mean_diploids,mean_triploids,mean_tetraploids,diploids,triploids,tetraploids,population_per_generation_per_rep

def fitness(individual1,individual2,selection_coeff=0.2,reference_fitness=0.8):
    """Implemented fitness function (see manuscript)

    Args:
        individual1 (lst): list containing genetic information which forms the individual
        individual2 (_type_): list containing genetic information which forms the individual
        selection_coeff (float, optional): selection coefficient used in formula. Defaults to 0.2.
        reference_fitness (float, optional):  reference fitness used in formula. Defaults to 0.8.

    Returns:
        Boolean value that determines if mating proceeds or not
    """
    p1, p2 = len(individual1), len(individual2)
    # Allele A increases fitness and allele B decreases fitness
    All_alleles_ind_1,All_alleles_ind_2 = [element for sublist in individual1 for element in sublist], [element for sublist in individual2 for element in sublist]
    fitness_score1 = reference_fitness + ((All_alleles_ind_1.count("A")/p1) * selection_coeff) - ((All_alleles_ind_1.count("B")/p1) * selection_coeff)
    fitness_score2 = reference_fitness + ((All_alleles_ind_2.count("A")/p2) * selection_coeff) - ((All_alleles_ind_2.count("B")/p2) * selection_coeff)
    # check if mating fails
    random_n1,random_n2 = random.uniform(0,1),random.uniform(0,1)
    scores = [fitness_score1,fitness_score2]
    random_n = [random_n1,random_n2]
    if max(random_n) < min(scores):
        return True
    else:
        return False

def allelefreq_per_generation_average_cytotypespecific(population_per_generation_per_rep, ploidy_level, freq=True):
    """Function that plots the allelefrequencies in each generation per cytotype

    Args:
        population_per_generation_per_rep (lst): lst containing the population (all individuals) of each generation per repetition
        ploidy_level (int): number of chromosome sets
        freq (bool, optional): _description_. Defaults to True.
    """
    average_A = []
    average_a = []
    average_B = []
    average_b = []
    for rep in population_per_generation_per_rep:  # Loop over each individual repetition
        freq_A = []
        freq_a = []
        freq_B = []
        freq_b = []
        for gen_pop in rep: # loop over every generation in each repetition
            count_A = []
            count_a = []
            count_B = []
            count_b = []
            cytotype_specific = []
            for individual in gen_pop:
                if len(individual) == ploidy_level: # select only the cytotype of interest
                    cytotype_specific.append(individual)
            for gen in cytotype_specific:                     # count the number of each specific allele present per generation in the specific cytotype
                gen = [element for sublist in gen for element in sublist]
                gen = [element for sublist in gen for element in sublist]
                count_A.append(gen.count("A"))
                count_a.append(gen.count("a"))
                count_B.append(gen.count("B"))
                count_b.append(gen.count("b"))
            if freq == True:
                amount_of_alleles = int(sum(count_A) + sum(count_a) + sum(count_B) + sum(count_b))
                if amount_of_alleles != 0: 
                    freq_A.append(sum(count_A)/amount_of_alleles)
                    freq_a.append(sum(count_a)/amount_of_alleles)
                    freq_B.append(sum(count_B)/amount_of_alleles)
                    freq_b.append(sum(count_b)/amount_of_alleles)
                else:
                    freq_A.append(0)
                    freq_a.append(0)
                    freq_B.append(0)
                    freq_b.append(0)
            else:
                freq_A.append(sum(count_A))
                freq_a.append(sum(count_a))
                freq_B.append(sum(count_B))
                freq_b.append(sum(count_b))
        average_A.append(freq_A)
        average_a.append(freq_a)
        average_B.append(freq_B)
        average_b.append(freq_b)
    average_A = get_mean_of_freq_overreps_counts(average_A,len(population_per_generation_per_rep))
    average_a = get_mean_of_freq_overreps_counts(average_a,len(population_per_generation_per_rep))
    average_B = get_mean_of_freq_overreps_counts(average_B,len(population_per_generation_per_rep))
    average_b = get_mean_of_freq_overreps_counts(average_b,len(population_per_generation_per_rep))
    generations = list(range(1,len(population_per_generation_per_rep[0]) + 1))
    fig= plt.figure()
    fig.suptitle("Allele frequencies per cytotype")
    sns.lineplot(x=generations, y=average_A, color="red")
    sns.lineplot(x=generations, y=average_a, color="blue")
    sns.lineplot(x=generations, y=average_B, color="green")
    sns.lineplot(x=generations, y=average_b, color="yellow")
    if freq == True:
        plt.ylim(0,1)
    plt.show()

def fitness_of_population(population_per_generations_per_rep, reference_fitness, selection_coeff):
    """Function that determines the average fitness per population per generation

    Args:
        population_per_generation_per_rep (lst): lst containing the population (all individuals) of each generation per repetition
        reference_fitness (float): selection coefficient
        selection_coeff (float): selection coefficient

    Returns:
        returns a list of fitness score of each population per generation per repetition
    """
    scores_per_rep = []
    for rep in population_per_generations_per_rep: # loop over each repetition
        score_per_generation = [] 
        for pop_per_generation in rep: # loop over all generations which contain the population in each generation
            score_per_generation.append(fitnesscore_per_population(pop_per_generation,reference_fitness=reference_fitness,selection_coeff=selection_coeff))
        scores_per_rep.append(score_per_generation)
    mean_fitness_of_population_over_reps = get_mean_of_freq_overreps_counts(scores_per_rep,len(scores_per_rep))
    return mean_fitness_of_population_over_reps

def fitnesscore_per_population(population,reference_fitness=0.8,selection_coeff=0.2):
    """Function that determines the fitnesss score of a population

    Args:
        population (lst): list of individuals
        reference_fitness (float, optional): reference fitness score Defaults to 0.8.
        selection_coeff (float, optional): selection coefficient. Defaults to 0.2.

    Returns:
        Fitness score of the entire population
    """
    score = []
    for individual in population:
        All_alleles_ind = [element for sublist in individual for element in sublist]
        fitness_score = reference_fitness + ((All_alleles_ind.count("A")/len(individual)) * selection_coeff) - ((All_alleles_ind.count("B")/len(individual)) * selection_coeff)
        score.append(fitness_score)
    return sum(score)/len(population)

def plot_fitness_vs_gen(fitness_per_generation):
    """Function to plot the fitness of the population in each generation

    Args:
        fitness_per_generation (lst): list containing the fitness score of the population in each generation
    """
    generations = list(range(1,len(fitness_per_generation)+1))
    fig= plt.figure()
    fig.suptitle("Population fitness per generation")
    sns.lineplot(x=generations, y=fitness_per_generation, color="red")
    #plt.ylim(0.8,1)
    plt.show()

def cytotype_dynamics_buffer_experiment(initial_size:int,ploidy_number:int,loci_number:int,max_generations:int,v:float,f:float,reps:int,apply_fitness=True,selection_coeff=0.2,reference_fitness=0.8, recombination_frequency=0.01,type=None,amount_of_mutations=1,experiment=0):
    """The individual based stochastic simulation used in the buffer experiment as described in the manuscript

    Args:
        initial_size (int): initial population size (N)
        ploidy_number (int): ploidy level of initial population
        loci_number (int): amount of loci of each chromosome
        max_generations (int): maximum amount of generations the IBM should be run
        v (float): probability of unreduced gamete production
        f (float): triploid fitness
        reps (int): amount of repetitions the IBM should be run for
        apply_fitness (bool, optional): apply fitness function or not. Defaults to True.
        selection_coeff (float, optional): selection coefficient of fitness. Defaults to 0.2.
        reference_fitness (float, optional): reference fitness score when fitness is applied. Defaults to 0.8.
        recombination_frequency (float, optional): probability that recombination takes place per loci. Defaults to 0.01.
        type (_type_, optional): type of population. Defaults to None.
        amount_of_mutations (int, optional): amount of individuals introduced with beneficial mutation. Defaults to 1.
        experiment (int, optional): to see if buffer experiment happens. Defaults to 0.

    Returns:
        amount of generations it takes for the population to reach the maximum fitness of 1
    """

    results = []
    for i in range(0,reps):
        initial_population = initialize_population(initial_size,loci_number,ploidy_number,type=type) # na each rep opnieuw initialpop opstarten ==> dus een nieuwe onafhankelijke trials
        t = 1
        n = len(initial_population)
        new_pop = initial_population
        population_per_generation = [initial_population]
        while t < max_generations:
            if experiment == 0:
                if fitnesscore_per_population(new_pop) == 1:
                    results.append([t-1,amount_of_mutations])
                    break
            else:
                results.append(fitnesscore_per_population(new_pop))
            #print(fitnesscore_per_population(new_pop))
            offsprings = []
            while len(offsprings) < n:
                if t == 150:
                    if len(offsprings) == 0:
                    # add individual with mutation
                        a = 1
                        while a < amount_of_mutations:
                            if ploidy_number == 2:
                                offsprings.append([["A","b"],["A","b"]])
                            if ploidy_number == 4:
                                offsprings.append([["A","b"],["A","b"],["A","b"],["A","b"]])
                            a += 1
                pair = random.sample(new_pop,2)
                if apply_fitness:
                    if fitness(pair[0],pair[1],selection_coeff=0.2,reference_fitness=0.8) == True:
                        first = meiosis(pair[0],v,f,recombination_frequency)
                        second = meiosis(pair[1],v,f,recombination_frequency)
                        if first != () and second != () and len(list(first + second)) < 5:
                            offspring = list(first + second)
                            offsprings.append(offspring)
                else:
                    first = meiosis(pair[0],v,f,recombination_frequency)
                    second = meiosis(pair[1],v,f,recombination_frequency)
                    if first != () and second != () and len(list(first + second)) < 5:
                        offspring = list(first + second)
                        offsprings.append(offspring)                  
            new_pop = offsprings
            population_per_generation.append(offsprings)
            t += 1
    return results

def make_df_from_data(data,ploidy):
    mutations = []
    generations = []
    for point in data:
        if point != None:
            mutations.append(int(point[1]))
            generations.append(int(point[0]))
    name_ploidy = [f"{ploidy}"] * len(mutations)
    structure = {"ploidy":name_ploidy,"mutations":mutations, "generations":generations}
    df = pd.DataFrame(structure)
    return df

def buffer_experiment():
    """Function to run and plot the buffer experiment

    Returns:
        results of experiment
    """
    # first for diploid starting population with "normal fitness" of 0.8
    # to eliminate formation of triploids use fertility factor of 0
    # for 3 kind of v: 0%, 5% and 10%
    results_di = [] # contain list of list for each v with mutations and generation that fitness 1 is reached
    for v in [0.1]:
        for i in [5,50,100,200,400,500]:
            results_di.append(cytotype_dynamics_buffer_experiment(1000,4,2,500,v,0.05,50,True,0.2,0.8,0.01,0,i,0))
    results = make_df_from_data([list for list in results_di for list in list],"diploid")
    grouped = results.groupby(['ploidy', 'mutations'])['generations']
    mean = grouped.mean()
    std = grouped.std()
    print(mean)
    print(std)
    return results

def plot_buffer_experiment(data):
    """Function to plot the buffer experiment

    Args:
        data (lst): maximum amount of generations necessary to reach maximum polulation fitness for each number of mutations added
    """
    sns.boxplot(data=data, x="mutations", y="generations", hue="ploidy")
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
    #plt.savefig('bufferexperiment.png')
    plt.show()



# figure 5
#phase_portrait() 

# figure 6
#plot_fertility_unreduced_gametes()

#Figure 7
#cytotypespecific_frequencies_di_tri_tetraploid(cytotype_dynamics(1000,2,2,max_generations=500,v=0.01,f=0.3,reps=10,apply_fitness=False,type=None),recursive_NL_equations(max_generations=500,v=0.01,f=0.3,d1=0.25,d2=0.25,d3=0.5))

# Figure 8 and 9
#large_pop = cytotype_dynamics(1000,2,2,500,0.01,0.3,30,True)
#cytotypespecific_frequencies_di_tri_tetraploid(large_pop,recursive_NL_equations(500,0.01,0.3,0.25,0.25,0.5))

# Figure 10 Check allele count and frequencies in each cytotype
#large_pop = cytotype_dynamics(1000,2,2,500,0,0.3,30,True)
#allelefreq_per_generation_average_cytotypespecific(large_pop[-1],3,True)

# Figure 11
#plot_fitness_vs_gen(fitness_of_population(cytotype_dynamics(1000,2,2,500,0,0.3,30,True,recombination_frequency=0.1)[-1],0.8,0.2))

# Figure 11-15
# plot_buffer_experiment(buffer_experiment())