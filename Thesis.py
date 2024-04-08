# Masterthesis
import random
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

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

def cytotype_dynamics(initial_size:int,ploidy_number:int,loci_number:int,max_generations:int,v:float,f:float,reps:int,apply_fitness=True,selection_coeff=0.2,reference_fitness=0.8, recombination_frequency=0.01):
    diploids = []
    triploids = []
    tetraploids = []
    population_per_generation_per_rep = []
    for i in range(0,reps):
        initial_population = initialize_population(initial_size,loci_number,ploidy_number) # na each rep opnieuw initialpop opstarten ==> dus een nieuwe onafhankelijke trials
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
    return max_generations,mean_diploids,mean_triploids,mean_tetraploids,diploids,triploids,tetraploids,population_per_generation_per_rep[0]

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
    if min(random_n) < min(scores):
        return False
    else:
        return True

def IBM_cytotypespecific_frequencies_di_tri_tetraploid(data_simulation):
    t = list(range(1,data_simulation[0]+1))
    fig,axes = plt.subplots(3,1,sharex=True)
    fig.suptitle("Cytotype-specific dynamics")
    axes.set(xlabel='Generations', ylabel='Frequency')
    # Diploid
    for rep in data_simulation[4]:
        sns.scatterplot(ax= axes[0],x=t, y=rep, color="red")
    plt.ylim(0, 1)
    # Triploid
    for rep in data_simulation[5]:
        sns.scatterplot(ax= axes[1], x=t, y=rep, color="green")
    plt.ylim(0, 0.8)
    # Tetraploid 
    for rep in data_simulation[6]:
        sns.scatterplot(ax= axes[2], x=t, y=rep, color="blue")
    plt.ylim(0, 0.1)
    plt.show()

large_pop = cytotype_dynamics(100,2,2,500,0.01,0.3,1,True)
print(len(large_pop[-1]))

def allelefreq_per_generation(population_per_generation):
    # divide per cytotype
    diploids = []
    triploids = []
    tetraploids = []
    for pop in population_per_generation:
        for individual in pop:
            if len(individual) == 2:
                diploids.append(individual)
            if len(individual) == 3:
                triploids.append(individual)
            if len(individual) == 4:
                tetraploids.append(individual)
    return diploids,tetraploids,tetraploids

print(len(allelefreq_per_generation(large_pop[-1])[0]))
            






# Idea: Plot frequencies of each allele per generation per for each cytotype



# figure 2
#plot_fertility_unreduced_gametes()

# figure 3
#phase_portrait() 
    
#Figure 4
#cytotypespecific_frequencies_di_tri_tetraploid(cytotype_dynamics(500,2,2,max_generations=500,v=0.01,f=0,reps=10),recursive_NL_equations(max_generations=500,v=0.01,f=0,d1=0.25,d2=0.25,d3=0.5))

# Figure 4 and 5
#cytotypespecific_frequencies_di_tri_tetraploid(large_pop,recursive_NL_equations(500,0.01,0.3,0.25,0.25,0.5))