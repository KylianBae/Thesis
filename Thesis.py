# Masterthesis
import random
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

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

def plot_line(inp,title):
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
    # haploid gamete frequency: 2g^2(2-v-2f) + g(-3+3v+3f) + (1-v-f) = 0         
    # change in frequency: x * (4 - 2 * v - 4 * f) + (-3 + 3 * v1 + 3 * f)           
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
        5= population of triploids with worst genotype


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
    mean_values = [0] * len(data[0])
    for rep in data:
        for i, value in enumerate(rep):
            mean_values[i] += value
    mean_values = [(value/reps) for value in mean_values]
    return mean_values

def get_mean_of_freq_overreps_counts(data,reps):
    mean_values = [0] * len(data[0])
    for rep in data:
        for i, value in enumerate(rep):
            if not isinstance(value,list):
                mean_values[i] += value
    mean_values = [(value/reps) for value in mean_values]
    return mean_values

def recombination(individual, recombination_frequency = 0.01, length = 2):
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
    scores_per_rep = []
    for rep in population_per_generations_per_rep: # loop over each repetition
        score_per_generation = [] 
        for pop_per_generation in rep: # loop over all generations which contain the population in each generation
            score_per_generation.append(fitnesscore_per_population(pop_per_generation,reference_fitness=reference_fitness,selection_coeff=selection_coeff))
        scores_per_rep.append(score_per_generation)
    mean_fitness_of_population_over_reps = get_mean_of_freq_overreps_counts(scores_per_rep,len(scores_per_rep))
    return mean_fitness_of_population_over_reps

def fitnesscore_per_population(population,reference_fitness=0.8,selection_coeff=0.2):
    score = []
    for individual in population:
        All_alleles_ind = [element for sublist in individual for element in sublist]
        fitness_score = reference_fitness + ((All_alleles_ind.count("A")/len(individual)) * selection_coeff) - ((All_alleles_ind.count("B")/len(individual)) * selection_coeff)
        score.append(fitness_score)
    return sum(score)/len(population)

def plot_fitness_vs_gen(fitness_per_generation):
    generations = list(range(1,len(fitness_per_generation)+1))
    fig= plt.figure()
    fig.suptitle("Population fitness per generation")
    sns.lineplot(x=generations, y=fitness_per_generation, color="red")
    #plt.ylim(0.8,1)
    plt.show()

# figure 2
#plot_fertility_unreduced_gametes()

# figure 3
#phase_portrait() 
    
#Figure 4
#cytotypespecific_frequencies_di_tri_tetraploid(cytotype_dynamics(1000,2,2,max_generations=500,v=0.01,f=0.3,reps=10,apply_fitness=False,type=None),recursive_NL_equations(max_generations=500,v=0.01,f=0.3,d1=0.25,d2=0.25,d3=0.5))

# Figure 4 and 5
#large_pop = cytotype_dynamics(1000,2,2,500,0.01,0.3,30,True)
#cytotypespecific_frequencies_di_tri_tetraploid(large_pop,recursive_NL_equations(500,0.01,0.3,0.25,0.25,0.5))

# Figure 7 Check allele count and frequencies in each cytotype
#large_pop = cytotype_dynamics(1000,2,2,500,0,0.3,30,True)
#allelefreq_per_generation_average_cytotypespecific(large_pop[-1],3,True)

# Figure 8
#plot_fitness_vs_gen(fitness_of_population(cytotype_dynamics(1000,2,2,500,0,0.3,30,True,recombination_frequency=0.1)[-1],0.8,0.2))

# Figure 9 effect fo cytotype buffering of mutations?
def cytotype_dynamics_buffer_experiment(initial_size:int,ploidy_number:int,loci_number:int,max_generations:int,v:float,f:float,reps:int,apply_fitness=True,selection_coeff=0.2,reference_fitness=0.8, recombination_frequency=0.01,type=None,amount_of_mutations=1):
    for i in range(0,reps):
        initial_population = initialize_population(initial_size,loci_number,ploidy_number,type=type) # na each rep opnieuw initialpop opstarten ==> dus een nieuwe onafhankelijke trials
        t = 1
        n = len(initial_population)
        new_pop = initial_population
        population_per_generation = [initial_population]
        while t < max_generations:
            if fitnesscore_per_population(new_pop) == 1:
                return [t-1,amount_of_mutations]
            #print(fitnesscore_per_population(new_pop))
            offsprings = []
            while len(offsprings) < n:
                if t == 100:
                    if len(offsprings) == 0:
                    # add individual with mutation
                        a = 1
                        while a < amount_of_mutations:
                            offsprings.append([["A","b"],["A","b"]])
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

    

#pop = cytotype_dynamics_buffer_experiment(1000,2,2,500,0.01,0.3,1,True,0.2,0.8,0.01,type=0)
# In buffered experiment of cytotpe dynamics we start with a population consisting of purely the worst possible genotype in terms of 
# fitness, we look at when the fitness of the dynamics simulation stabilizes and when it is stabilized we introduced an individual with a
# mutation which increases the fitness aka an A allele
#plot_fitness_vs_gen(fitness_of_population(pop[-1],0.8,0.2))

def buffer_experiment():
    di_generation_nmutations_coordinates = []
    for i in range(1,1001):
        di_generation_nmutations_coordinates.append(cytotype_dynamics_buffer_experiment(1000,2,2,500,0.01,0.3,1,True,0.2,0.8,0.01,0,i))
        print(f"Mutations: {i}")
    tri_generation_nmutations_coordinates = []
    for i in range(1,1001):
        tri_generation_nmutations_coordinates.append(cytotype_dynamics_buffer_experiment(1000,2,2,500,0.01,0.3,1,True,0.2,0.8,0.01,5,i))
        print(f"Mutations: {i}")
    tetra_generation_nmutations_coordinates = []
    for i in range(1,1001):
        tetra_generation_nmutations_coordinates.append(cytotype_dynamics_buffer_experiment(1000,2,2,500,0.01,0.3,1,True,0.2,0.8,0.01,4,i))
        print(f"Mutations: {i}")
    return di_generation_nmutations_coordinates,tri_generation_nmutations_coordinates,tetra_generation_nmutations_coordinates


def make_df_from_data(data_diploid, data_triploid, data_tetraploid):
    mutations_diploid = []
    generations_diploid = []
    for point in data_diploid:
        if point != None:
            mutations_diploid.append(int(point[1]))
            generations_diploid.append(int(point[0]))
    name_ploidy_di = ["Diploid"] * len(mutations_diploid)
    structure_diploid = {'ploidy': name_ploidy_di, "mutations":mutations_diploid, "generations":generations_diploid}
    df_diploid = pd.DataFrame(structure_diploid)

    mutations_triploid = []
    generations_triploid = []
    for point in data_triploid:
        if point != None:
            mutations_triploid.append(int(point[1]))
            generations_triploid.append(int(point[0]))
    name_ploidy_tri = ["Triploid"] * len(mutations_triploid)
    structure_triploid = {'ploidy': name_ploidy_tri, "mutations":mutations_triploid, "generations":generations_triploid}
    df_triploid = pd.DataFrame(structure_triploid)

    mutations_tetraploid = []
    generations_tetraploid = []
    for point in data_tetraploid:
        if point != None:
            mutations_tetraploid.append(int(point[1]))
            generations_tetraploid.append(int(point[0]))
    name_ploidy_tetra = ["Tetraploid"] * len(mutations_tetraploid)
    structure_tetra = {'ploidy': name_ploidy_tetra, "mutations":mutations_tetraploid, "generations":generations_tetraploid}
    df_tetra = pd.DataFrame(structure_tetra)
    
    final=pd.concat([df_diploid,df_triploid,df_tetra])
    return final

def plot_buffer_experiment(data_diploid,data_triploid,data_tetraploid):
    # split the data in mutations and corresponding generation that fitness 1 is achieved
    data = make_df_from_data(data_diploid,data_triploid,data_tetraploid)
    sns.jointplot(data=data, x="generations", y="mutations",hue="ploidy", s=10, linewidth=0)
    plt.show()


#data = [None, None, None, None, None, None, None, None, [429, 9], None, None, None, None, None, [381,15],[264,16],None, None, [355, 19], [427, 20], [377, 21], [332, 22], [450, 23], [462, 24], [419, 25], [435, 26], [249, 27], [335, 28], [419, 29], [294, 30], [492, 31], None, None, [429, 34], [365, 35], [372, 36], [296, 37], [257, 38], [297, 39], [336, 40], None, [344, 42], [492, 
#43], [430, 44], None, [323, 46], [407, 47], [312, 48], [416, 49], [244, 50], [200, 51], [431, 52], [264, 53], [471, 54], [390, 55], [446, 56], [461, 57], [471, 58], None, None, [413, 61], [376, 62], [288, 63], [364, 64], [400, 65], None, [431, 67], [371, 68], [443, 69], [272, 70], None, [344, 72], [409, 73], [298, 74], [315, 75], [190, 76], [436, 77], [181, 78], [205, 79], [288, 80], [342, 81], [238, 82], [209, 83], [444, 84], [259, 85], [195, 86], [315, 87], [358, 88], [386, 89], [454, 90], [396, 91], [219, 92], [439, 93], [318, 94], [406, 95], [313, 96], None, [278, 98], [293, 99], [210, 100], [213, 101], [310, 102], [262, 103], [195, 104], None, [251, 106], [335, 107], [314, 108], [285, 109], [436, 110], [163, 111], [202, 112], [174, 113], [183, 114], [274, 115], [243, 116], [183, 117], [292, 118], [268, 119], [219, 120], [170, 121], [281, 122], [160, 123], [176, 124], [177, 125], [171, 126], [300, 127], [285, 128], [435, 129], [210, 130], [211, 131], [204, 132], [302, 133], [178, 134], [394, 135], [260, 136], [315, 137], [207, 138], [236, 139], [178, 140], [189, 141], [195, 142], [262, 143], [249, 144], [176, 145], [257, 146], [286, 147], [188, 148], [144, 149], [192, 150], [201, 151], [352, 152], [265, 153], [223, 154], [157, 155], [357, 156], [305, 157], [159, 158], [164, 159], [195, 160], [197, 161], [151, 162], [141, 163], [154, 164], [191, 165], [157, 166], [198, 167], [158, 168], [218, 169], [148, 170], [175, 171], [208, 172], [221, 173], [227, 174], [232, 175], [377, 176], [190, 177], [328, 178], [183, 179], [483, 180], [200, 181], [204, 182], [154, 183], [188, 184], [172, 185], [173, 186], [166, 187], [152, 188], [179, 189], [142, 190], [207, 191], [215, 192], [144, 193], [198, 194], [225, 195], [211, 196], [238, 197], [165, 198], [183, 199], [151, 200], [184, 201], [179, 202], [201, 203], [165, 204], [196, 205], [140, 206], [169, 207], [209, 208], [216, 209], [182, 210], [198, 211], [173, 212], [173, 213], [223, 214], [274, 215], [168, 216], [153, 217], [155, 218], [348, 219], [166, 220], [188, 221], [145, 222], [203, 223], [163, 224], [158, 225], [252, 226], [384, 227], [184, 228], [198, 229], [183, 230], [171, 231], [167, 232], [211, 233], [151, 234], [149, 235], [290, 236], [158, 237], [198, 238], [163, 239], [155, 240], [165, 241], [143, 242], [178, 243], [162, 244], [194, 245], [143, 246], [168, 247], [160, 248], [223, 249], [142, 250], [268, 251], [191, 252], [211, 253], [182, 254], [397, 255], [149, 256], [162, 257], [157, 258], [212, 259], [230, 260], [202, 261], [143, 262], [155, 263], [166, 264], [303, 265], [163, 266], [144, 267], [249, 268], [238, 269], [149, 270], [186, 271], [149, 272], [201, 273], [160, 274], [137, 275], [148, 276], [155, 277], [182, 278], [150, 279], [168, 280], [228, 281], [144, 282], [195, 283], [135, 284], [162, 285], [143, 286], [151, 287], [159, 288], [150, 289], [155, 290], [148, 291], [189, 292], [177, 293], [274, 294], [158, 295], [173, 296], [162, 297], [225, 298], [140, 299], [160, 300], [289, 301], [161, 302], [188, 303], [146, 304], [124, 305], [159, 306], [161, 307], [256, 308], [251, 309], [185, 310], [166, 311], [161, 312], [147, 313], [184, 314], [158, 315], [141, 316], [164, 317], [164, 318], [166, 319], [120, 320], [200, 321], [173, 322], [145, 323], [190, 324], [180, 325], [153, 326], [163, 327], [144, 328], [168, 329], [171, 330], [163, 331], [170, 332], [233, 333], [154, 334], [209, 335], [146, 336], [148, 337], [165, 338], [146, 339], [134, 340], [139, 341], [145, 342], [178, 343], [137, 344], [225, 345], [210, 346], [142, 347], [258, 348], [144, 349], [166, 350], [160, 351], [164, 352], [137, 353], [170, 354], [204, 355], [163, 356], [165, 357], [125, 358], [174, 359], [127, 360], [145, 361], [142, 362], [133, 363], [149, 364], [147, 365], [171, 366], [129, 367], [278, 368], [144, 369], [156, 370], None, [146, 372], [131, 373], [307, 374], [187, 375], [169, 376], [153, 377], [194, 378], [155, 379], [229, 380], [175, 381], [202, 382], [159, 383], [133, 384], [244, 385], [144, 386], [191, 387], [182, 388], [178, 389], [163, 390], [173, 391], [156, 392], [156, 393], [150, 394], [193, 395], [134, 396], [152, 397], [148, 398], [150, 399], [190, 400], [279, 401], [190, 402], [171, 403], [140, 404], [152, 405], [146, 406], [142, 407], [162, 408], [203, 409], [156, 410], [133, 411], [260, 412], [140, 413], [222, 414], [154, 415], [146, 416], [148, 417], [152, 418], [152, 419], [154, 420], [153, 421], [160, 422], [138, 423], [135, 424], [241, 425], [142, 426], [133, 427], [141, 428], [134, 429], [151, 430], [155, 431], [150, 432], [158, 433], [132, 434], [180, 435], [142, 436], [132, 437], [132, 438], [130, 439], [170, 440], [168, 441], [145, 442], [147, 443], [150, 444], [128, 445], [165, 446], [152, 447], [135, 448], [168, 449], [159, 450], [134, 451], [137, 452], [167, 453], [183, 454], [179, 455], [156, 456], [139, 457], [151, 458], [140, 459], [128, 460], [199, 461], [191, 462], [134, 463], [135, 464], [143, 465], [191, 466], [191, 467], [142, 468], [141, 469], [304, 470], [144, 471], [194, 472], [128, 473], [140, 474], [145, 475], [171, 476], [207, 477], [143, 478], [130, 479], [120, 480], [167, 481], [159, 482], [197, 483], [158, 484], [133, 485], [241, 486], [157, 487], [139, 488], [172, 489], [159, 490], [144, 491], [159, 492], [144, 493], [126, 494], [179, 495], [122, 496], [259, 497], [184, 498], [168, 499], [137, 500], [139, 501], [135, 502], [133, 503], [124, 504], [161, 505], [133, 506], [139, 507], [123, 508], [127, 509], [124, 510], [142, 511], [139, 512], [142, 513], [192, 514], [142, 515], [130, 516], [128, 517], [260, 518], [214, 519], [127, 520], [160, 521], [169, 522], [128, 523], [154, 524], [133, 525], [136, 526], [226, 527], [151, 528], [131, 529], [132, 530], [211, 531], [173, 532], [133, 533], [123, 534], [133, 535], [139, 536], [134, 537], [131, 538], [156, 539], [125, 540], [137, 541], [320, 542], [166, 543], [163, 544], [121, 545], [140, 546], [158, 547], [139, 548], [151, 549], [284, 550], [163, 551], [122, 552], [188, 553], [290, 554], [133, 555], [160, 556], [126, 557], [127, 558], [230, 559], [156, 560], [150, 561], [155, 562], [145, 563], [164, 564], [144, 565], [154, 566], [152, 567], [140, 568], [142, 569], [123, 570], [130, 571], [142, 572], [132, 573], [145, 574], [148, 575], [135, 576], [132, 577], [197, 578], [155, 579], [147, 580], [140, 581], [144, 582], [169, 583], [158, 584], [141, 585], [122, 586], [147, 587], [140, 588], [161, 589], [144, 590], [169, 591], [153, 592], [146, 593], [177, 594], [152, 595], [149, 596], [161, 597], [138, 598], [141, 599], [172, 600], [144, 601], [112, 602], [146, 603], [137, 604], [134, 605], [129, 606], [129, 607], [123, 608], [122, 609], [120, 610], [121, 611], [123, 612], [130, 613], [122, 614], [137, 615], [149, 616], [148, 617], [140, 618], [125, 619], [130, 620], [142, 621], [159, 622], [143, 623], [132, 624], [139, 625], [136, 626], [188, 627], [144, 628], [132, 629], [133, 630], [151, 631], [135, 632], [129, 633], [147, 634], [127, 635], [135, 636], [157, 637], [126, 638], [152, 639], [140, 640], [123, 641], [149, 642], [119, 643], [123, 644], [121, 645], [138, 646], [127, 647], [127, 648], [138, 649], [128, 650], [151, 651], [136, 652], [143, 653], [136, 654], [146, 655], [127, 656], [144, 657], [130, 658], [120, 659], [129, 660], [147, 661], [127, 662], [161, 663], [120, 664], [140, 665], [119, 666], [148, 667], [134, 668], [148, 669], [140, 670], [121, 671], [130, 672], [126, 673], [120, 674], [126, 675], [121, 676], [138, 677], [123, 678], [161, 679], [203, 680], [137, 681], [125, 682], [123, 683], [135, 684], [124, 685], [115, 686], [135, 687], [134, 688], [136, 689], [118, 690], [133, 691], [121, 692], [122, 693], [129, 694], [147, 695], [125, 696], [118, 697], [120, 698], [123, 699], [132, 700], [168, 701], [153, 702], [146, 703], [132, 704], [131, 705], [135, 706], [121, 707], [177, 708], [121, 709], [132, 710], [176, 711], [135, 712], [129, 713], [126, 714], [129, 715], [120, 716], [131, 717], [128, 718], [144, 719], [123, 720], [168, 721], [134, 722], [128, 723], [240, 724], [141, 725], [164, 726], [125, 727], [141, 728], [157, 729], [124, 730], [125, 731], [134, 732], [125, 733], [122, 734], [124, 735], [119, 736], [161, 737], [130, 738], [126, 739], [132, 740], [114, 741], [126, 742], [113, 743], [128, 744], [145, 745], [115, 746], [126, 747], [120, 748], [133, 749], [115, 750], [127, 751], [141, 752], [184, 753], [141, 754], [125, 755], [130, 756], [144, 757], [257, 758], [245, 759], [136, 760], [129, 761], [131, 762], [113, 763], [140, 764], [127, 765], [116, 766], [167, 767], [120, 768], [124, 769], [122, 770], [166, 771], [128, 772], [138, 773], [118, 774], [112, 775], [124, 776], [119, 777], [125, 778], [160, 779], [124, 780], [140, 781], [128, 782], [174, 783], [127, 784], [130, 785], [125, 786], [132, 787], [126, 788], [138, 789], [116, 790], [127, 791], [126, 792], [137, 793], [114, 794], [120, 795], [154, 796], [120, 797], [143, 798], [159, 799], [137, 800], [140, 801], [119, 802], [127, 803], [114, 804], [116, 805], [119, 806], [114, 807], [126, 808], [162, 809], [118, 810], [122, 811], [114, 812], [126, 813], [133, 814], [126, 815], [113, 816], [123, 817], [128, 818], [115, 819], [129, 820], [134, 821], [118, 822], [121, 823], [123, 824], [112, 825], [118, 826], [117, 827], [126, 828], [140, 829], [125, 830], [123, 831], [127, 832], [139, 833], [114, 834], [136, 835], [123, 836], [115, 837], [131, 838], [148, 839], [123, 840], [112, 841], [112, 842], [144, 843], [111, 844], [132, 845], [124, 846], [115, 847], [121, 848], [122, 849], [128, 850], [129, 851], [123, 852], [145, 853], [130, 854], [134, 855], [119, 856], [139, 857], [170, 858], [127, 859], [116, 860], [124, 861], [128, 862], [113, 863], [125, 864], [118, 865], [120, 866], [115, 867], [117, 868], [118, 869], [170, 870], [122, 871], [125, 872], [209, 873], [145, 874], [212, 875], [114, 876], [134, 877], [119, 878], [162, 879], [113, 880], [122, 881], [122, 882], [120, 883], [148, 884], [115, 885], [121, 886], [209, 887], [114, 888], [119, 889], [162, 890], [117, 891], [118, 892], [171, 893], [116, 894], [122, 895], [118, 896], [120, 897], [121, 898], [118, 899], [120, 900], [117, 901], [146, 902], [120, 903], [114, 904], [122, 905], [115, 906], [112, 907], [118, 908], [124, 909], [119, 910], [112, 911], [121, 912], [113, 913], [141, 914], [140, 915], [115, 916], [111, 917], [135, 918], [113, 919], [144, 920], [118, 921], [112, 922], [142, 923], [124, 924], [110, 925], [120, 926], [118, 927], [111, 928], [108, 929], [110, 930], [124, 931], [115, 932], [116, 933], [120, 934], [154, 935], [109, 936], [111, 937], [112, 938], [107, 939], [117, 940], [106, 941], [116, 942], [129, 943], [115, 944], [111, 945], [118, 946], [116, 947], [111, 948], [124, 949], [134, 950], [113, 951], [119, 952], [112, 953], [109, 954], [122, 955], [113, 956], [130, 957], [109, 958], [107, 959], [127, 960], [109, 961], [128, 962], [132, 963], [117, 964], [110, 965], [119, 966], [129, 967], [116, 968], [111, 969], [111, 970], [106, 971], [112, 972], [124, 973], [106, 974], [110, 975], [125, 976], [116, 977], [115, 978], [109, 979], [108, 980], [116, 981], [105, 982], [135, 983], [104, 984], [103, 985], [107, 986], [109, 987], [105, 988], [125, 989], [103, 990], [104, 991], [105, 992], [105, 993], [102, 994], [106, 995], [104, 996], [106, 997], [102, 998], [102, 999]]
data = buffer_experiment()
print(data)
plot_buffer_experiment(data[0],data[1],data[2])