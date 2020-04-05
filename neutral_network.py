import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from Levenshtein import distance as levenshtein_distance
from difflib import ndiff
import random

# Map from https://www.biostars.org/p/2903/
encodings = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
    "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
    "UAU":"Y", "UAC":"Y", "UAA":"STOP", "UAG":"STOP",
    "UGU":"C", "UGC":"C", "UGA":"STOP", "UGG":"W",
    "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
    "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
    "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
    "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",}

survivors = {"UACUACCUUUUCAACUAUCAACAAACAAAUUAC" : True}

covid_protein = "YYLFNYQQTNY"
sars_protein = "YYYLNYNYTTY"
bat_protein = "YYSFNYNYTNY"
civet_protein = "YYYLNYKYTSY"
covid_nucleotides = "UACUACCUUUUCAACUAUCAACAAACAAAUUAC"
draw_graph = True
show_failures = False

survival_probability = 0.15

def find_first_different_char(string1, string2):
    s1, s2 = list(string1), list(string2)

    for i in range(0, len(s1)):
        if s1[i] != s2[i]:
            return i

    return -1


def protein_from_nucleotides(genome):
    protein = ""
    for i in range(0,int(len(genome) / 3)):
        protein = protein + encodings[genome[3*i:(3*i)+3]]

    return protein

def get_key(val): 
    for key, value in encodings.items(): 
         if val == value: 
             return key 
  
    return "key doesn't exist"

def mutate(genome):
    new_genome = list(genome)
    point = random.Random().randrange(len(genome))
    bases = ['U', 'C', 'A', 'G']
    bases.remove(genome[point])
    base = bases[random.Random().randrange(0,3)]
    new_genome[point] = base
    return ''.join(new_genome)

def mutate_towards_target(genome, protein):
    new_genome = list(genome)
    protein_list = list(protein)
    current_protein = protein_from_nucleotides(genome)
    pro = find_first_different_char(current_protein, protein)
    ng0 = new_genome[3*pro:(3*pro)+3]
    first_mistake = find_first_different_char(ng0, get_key(protein_list[pro]))
    if pro == -1:
        point = random.Random().randrange(len(genome))
    else:
        point = random.Random().randrange(pro*3, pro*3 + 3)
    bases = ['U', 'C', 'A', 'G']
    bases.remove(genome[point])
    base = bases[random.Random().randrange(0,3)]
    new_genome[point] = base
    ng = new_genome[3*pro:(3*pro)+3]
    new_first_mistake = find_first_different_char(ng, get_key(protein_list[pro]))
    if new_first_mistake > first_mistake or new_first_mistake == -1:
        survivors[''.join(new_genome)] = True
    return ''.join(new_genome)

def determine_survival(genome):
    protein = protein_from_nucleotides(genome)
    for g in survivors:
        if protein_from_nucleotides(g) == protein:
            survivors[genome] = survivors[g]
            return survivors[g]
        
    if random.Random().random() < survival_probability:
        return True

    return False

network_labels = ["UACUACCUUUUCAACUAUCAACAAACAAAUUAC"]

colors = ['#00FF00']

if draw_graph:
    G = nx.Graph()

j = 0
danger_mutations = 0

for i in range(1,2000):
    if j == 0:
        orig = 0
    else:
        orig = random.Random().randrange(0, j)
        while not survivors[network_labels[orig]]:
            orig = random.Random().randrange(0, j)
    genome = mutate(network_labels[orig])
    #genome = mutate_towards_target(network_labels[orig], sars_protein)
    if genome in network_labels and draw_graph:
        if not show_failures:
            if genome in survivors and survivors[genome]:
                G.add_edge(orig, network_labels.index(genome))
        else:
            G.add_edge(orig, network_labels.index(genome))  
        continue
    j = j + 1
    network_labels.append(genome)
    if protein_from_nucleotides(genome) == bat_protein:
        if draw_graph:
            colors.append("#001177")
            G.add_edge(orig, j)
            survivors[genome] = True
            print("Became Bat")
    if protein_from_nucleotides(genome) == civet_protein:
        if draw_graph:
            colors.append("#3333CC")
            G.add_edge(orig, j)
            survivors[genome] = True
            print("Became Civet")
    if protein_from_nucleotides(genome) == sars_protein:
        if draw_graph:
            colors.append("#0000FF")
            G.add_edge(orig, j)
        print("Oh no! It became SARS")
        print(i)
        break
    if not network_labels[j] in survivors:
        survivors[network_labels[j]] = determine_survival(genome)
    if not survivors[network_labels[j]]:
        if show_failures:
            if draw_graph:
                colors.append('#FF0000')
                G.add_edge(orig,j)

    else:
        if draw_graph:
            if levenshtein_distance(protein_from_nucleotides(genome), covid_protein) > 5:
                colors.append('#FF00FF')
                danger_mutations = danger_mutations + 1
            else:
                colors.append('#AAAAAA')
            G.add_edge(orig, j)

while len(colors) < G.number_of_nodes():
    colors.insert(2, '#AAAAAA')

print(danger_mutations)
print(G.number_of_nodes())

if draw_graph:
    nx.draw(G, node_color = colors)
    plt.show()
