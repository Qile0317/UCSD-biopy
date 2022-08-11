from math import *
import pip
import pandas
import numpy as np 
import sys
import random

#There may be some problem with this function
def PatternCount(symbol,ExtendedGenome ):
    count = 0
    for i in range(0,len(ExtendedGenome)-len(symbol)+1):
        if ExtendedGenome[i:i+len(symbol)] == symbol:
            count = count+1
    return count   

#Testing function
Text = "CGCGATACGTTACATACATGATAGACCGCGCGCGATCATATCGCGATTATC"

Pattern = "CGCG" 

print(PatternCount(Text, Pattern))

#FrequencyMap for loop
def FrequencyMap(Text, k):
    freq = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        freq[Pattern] = 0
        for i in range(n-k+1):
            if Text[i:i+k] == Pattern:
                freq[Pattern] = freq[Pattern] + 1
    return freq

#Frequentwords
def FrequentWords(Text, k):
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())
    for key in freq:
        if freq[key] == m:
            p = key
            words.append(p)
    return words

# Input:  A string Pattern
# Output: The reverse of Pattern
def Reverse(Pattern):
    reverse = ''
    for char in Pattern:
        reverse = char + reverse
    return reverse

# Input:  A DNA string Pattern
# Output: The complementary string of Pattern (with every nucleotide replaced by its complement).
def Complement(Pattern):
    comp = " "
    basepair = {"A":"T", "T":"A", "C":"G", "G":"C"}
    for char in Pattern:
        comp += basepair.get(char)        #get() returns the value for the specified key
    return comp

def ReverseComplement(Pattern):
    return Reverse(Complement(Pattern))

#finding a pattern in string that matches
def PatternMatching(Pattern, Genome):
    positions = [] # output variable
    n = len(Genome)
    k = len(Pattern)
    for i in range(n-k+1):
        if Genome[i:i+k] == Pattern:
            positions.append(i)
    return positions

#Symbol array, dependent on patterncount
def SymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    for i in range(n):
        array[i] = PatternCount(symbol, ExtendedGenome[i:i+(n//2)])
    return array 

#more efficient version!
def FasterSymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    array[0] = PatternCount(symbol, Genome[0:n//2])
    for i in range(1, n):
        array[i] = array[i-1]
        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i]-1
        if ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1
    return array

#with open('e_coli.txt') as file:
    #e_coli = file.read();

#array = FasterSymbolArray(e_coli, "C")

#plt.plot(*zip(*sorted(array.items())))
#plt.show()

# Input:  A String Genome
# Output: The skew array of Genome as a list.
def SkewArray(Genome):
    skew = [0]
    score = {"A":0, "T":0, "C":-1, "G":1}
    for i in range(1,len(Genome)+1):
            skew.append(score[Genome[i-1]] + skew[i-1])
    return skew

#example with plot:
#stest = [0, -1, -1, -1, 0, 1, 2, 1, 1, 1, 0, 1, 2, 1, 0, 0, 0, 0, -1, 0, -1, -2]
#plt.plot(stest[:], marker='o')
#plt.show()

#minArray
def MinimumSkew(Genome):
    Positions = []
    Array = SkewArray(Genome)
    m = min(Array)
    for i in range(len(Array)):
        if Array[i] == m:
            Positions.append(i)
    return Positions

#maxArray
def MaximumSkew(Genome):
    Positions = []
    Array = SkewArray(Genome)
    m = max(Array)
    for i in range(len(Array)):
        if Array[i] == m:
            Positions.append(i)
    return Positions

# Input:  Two strings p and q
# Output: An integer value representing the Hamming Distance between p and q.
def HammingDistance(p, q):
    n_p = len(p)
    hd = 0
    for i in range(n_p):
        if p[i] != q[i]:
            hd = hd + 1
    return hd

#test = 43:
HammingDistance("CTTGAAGTGGACCTCTAGTTCCTCTACAAAGAACAGGTTGACCTGTCGCGAAG",
 "ATGCCTTACCTAGATGCAATGACGGACGTATTCCTTTTGCCTCAACGGCTCCT")

# Input:  Strings Pattern and Text along with an integer d
# Output: A list containing all starting positions where Pattern appears
# as a substring of Text with at most d mismatches
def ApproximatePatternMatching(Text, Pattern, d):
    positions = [] # initializing list of positions
    for i in range(len(Text)-len(Pattern)+1):
        if HammingDistance(Text[i:i+len(Pattern)], Pattern) <= d:
            positions.append(i)
    return positions

#slightly modified to include a count of how many:
def ApproximatePatternCount(Pattern, Text, d):
    count = 0 # initialize count variable
    positions = [] # initializing list of positions
    for i in range(len(Text)-len(Pattern)+1):
        if HammingDistance(Text[i:i+len(Pattern)], Pattern) <= d:
            positions.append(i)
            count += 1
    return count

                    ############################
##################### BRUTE FORCE MOTIF METHOD #####################################
                    ############################

# Input:  Motifs, a list of kmers (strings)
# Output: Count(Motifs), the count matrix of Motifs as a dictionary of lists
def Count(Motifs):
    count = {} 
    k = len(Motifs[0]) 
    for symbol in "ACGT":
        count[symbol] = [] # count matrix now has keys A, C, T, and G all with values of empty list
        for j in range(k):
            count[symbol].append(0) # count matrix now has keys A, C, G, and T all with values of a list of zeroes of length equal to the length of a kmer
    t = len(Motifs) # length of Motifs, a list of kmers (strings)
    for i in range(t): # for each kmer in Motifs
        for j in range(k): # for each element of the kmer
            symbol = Motifs[i][j] # assigns the key (symbol) to a nucleotide (ACGT) in Motifs
            count[symbol][j] += 1 # adds 1 to the position in the list assigned to the key
    return count

# Input:  A list of kmers Motifs
# Output: the profile matrix of Motifs, as a dictionary of lists.
def Profile(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {}
    profile = Count(Motifs)
    for symbol in "ACGT":
        for j in range(k):
            profile[symbol][j] = profile[symbol][j]/t
    return profile

# Input:  A set of kmers Motifs
# Output: A consensus string of Motifs.
def Consensus(Motifs):
    counts = Count(Motifs)
    consensus = []
    maxCounts = []
    loopCount = 0
    for symbol in counts:
        for i in range(len(counts[symbol])):
            if(loopCount == 0):
                consensus.append(symbol)
                maxCounts.append(counts[symbol][i])
            else:
                count = counts[symbol][i]
                if count > maxCounts[i]:
                    maxCounts[i] = count
                    consensus[i] = symbol
        loopCount += 1
    consensusString = ''
    for char in consensus:
        consensusString += char
    return consensusString

# Input:  A set of k-mers Motifs
# Output: The score of these k-mers.
def Score(Motifs):
    consensus = Consensus(Motifs)
    count = 0
    for motif in Motifs:
        for index in range(len(motif)):
            if motif[index] != consensus[index]:
                count += 1
    return count
                    #################################################################################
##################### GREEDY MOTIF (Faster but actually gets very inaccurate depending on datasets) ################################
                    #################################################################################

# Input:  String Text and profile matrix Profile
# Output: Pr(Text, Profile)
def Pr(Text, Profile):
    pro = 1
    for i in range(len(Text)):
        pro = pro*Profile[Text[i]][i]
    return pro

# The profile matrix assumes that the first row corresponds to A, the second corresponds to C,
# the third corresponds to G, and the fourth corresponds to T.
# You should represent the profile matrix as a dictionary, 
# whose keys are 'A', 'C', 'G', and 'T' and whose values are lists of floats
def ProfileMostProbableKmer(Text, k, Profile):
    max_probability = 0
    most_probable = Text[:k]
    n = len(Text)
    for i in range(n - k + 1):
        kmer = Text[i: i + k]
        probability = Pr(kmer, Profile)
        if probability > max_probability:
            max_probability = probability
            most_probable = kmer
    return most_probable

def GreedyMotifSearch(Dna,k,t): #this is the actual greedy algorithm, which has 6 dependencies of the previous code
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for m in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][m:m+k])
        for j in range(1, t):
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

#Application of the code we just wrote on 10 of 25 differentially expressed genes to find most probably motifs:
# Copy the ten strings occurring in the hyperlinked DosR dataset below.
Dna = [
"GCGCCCCGCCCGGACAGCCATGCGCTAACCCTGGCTTCGATGGCGCCGGCTCAGTTAGGGCCGGAAGTCCCCAATGTGGCAGACCTTTCGCCCCTGGCGGACGAATGACCCCAGTGGCCGGGACTTCAGGCCCTATCGGAGGGCTCCGGCGCGGTGGTCGGATTTGTCTGTGGAGGTTACACCCCAATCGCAAGGATGCATTATGACCAGCGAGCTGAGCCTGGTCGCCACTGGAAAGGGGAGCAACATC",
"CCGATCGGCATCACTATCGGTCCTGCGGCCGCCCATAGCGCTATATCCGGCTGGTGAAATCAATTGACAACCTTCGACTTTGAGGTGGCCTACGGCGAGGACAAGCCAGGCAAGCCAGCTGCCTCAACGCGCGCCAGTACGGGTCCATCGACCCGCGGCCCACGGGTCAAACGACCCTAGTGTTCGCTACGACGTGGTCGTACCTTCGGCAGCAGATCAGCAATAGCACCCCGACTCGAGGAGGATCCCG",
"ACCGTCGATGTGCCCGGTCGCGCCGCGTCCACCTCGGTCATCGACCCCACGATGAGGACGCCATCGGCCGCGACCAAGCCCCGTGAAACTCTGACGGCGTGCTGGCCGGGCTGCGGCACCTGATCACCTTAGGGCACTTGGGCCACCACAACGGGCCGCCGGTCTCGACAGTGGCCACCACCACACAGGTGACTTCCGGCGGGACGTAAGTCCCTAACGCGTCGTTCCGCACGCGGTTAGCTTTGCTGCC",
"GGGTCAGGTATATTTATCGCACACTTGGGCACATGACACACAAGCGCCAGAATCCCGGACCGAACCGAGCACCGTGGGTGGGCAGCCTCCATACAGCGATGACCTGATCGATCATCGGCCAGGGCGCCGGGCTTCCAACCGTGGCCGTCTCAGTACCCAGCCTCATTGACCCTTCGACGCATCCACTGCGCGTAAGTCGGCTCAACCCTTTCAAACCGCTGGATTACCGACCGCAGAAAGGGGGCAGGAC",
"GTAGGTCAAACCGGGTGTACATACCCGCTCAATCGCCCAGCACTTCGGGCAGATCACCGGGTTTCCCCGGTATCACCAATACTGCCACCAAACACAGCAGGCGGGAAGGGGCGAAAGTCCCTTATCCGACAATAAAACTTCGCTTGTTCGACGCCCGGTTCACCCGATATGCACGGCGCCCAGCCATTCGTGACCGACGTCCCCAGCCCCAAGGCCGAACGACCCTAGGAGCCACGAGCAATTCACAGCG",
"CCGCTGGCGACGCTGTTCGCCGGCAGCGTGCGTGACGACTTCGAGCTGCCCGACTACACCTGGTGACCACCGCCGACGGGCACCTCTCCGCCAGGTAGGCACGGTTTGTCGCCGGCAATGTGACCTTTGGGCGCGGTCTTGAGGACCTTCGGCCCCACCCACGAGGCCGCCGCCGGCCGATCGTATGACGTGCAATGTACGCCATAGGGTGCGTGTTACGGCGATTACCTGAAGGCGGCGGTGGTCCGGA",
"GGCCAACTGCACCGCGCTCTTGATGACATCGGTGGTCACCATGGTGTCCGGCATGATCAACCTCCGCTGTTCGATATCACCCCGATCTTTCTGAACGGCGGTTGGCAGACAACAGGGTCAATGGTCCCCAAGTGGATCACCGACGGGCGCGGACAAATGGCCCGCGCTTCGGGGACTTCTGTCCCTAGCCCTGGCCACGATGGGCTGGTCGGATCAAAGGCATCCGTTTCCATCGATTAGGAGGCATCAA",
"GTACATGTCCAGAGCGAGCCTCAGCTTCTGCGCAGCGACGGAAACTGCCACACTCAAAGCCTACTGGGCGCACGTGTGGCAACGAGTCGATCCACACGAAATGCCGCCGTTGGGCCGCGGACTAGCCGAATTTTCCGGGTGGTGACACAGCCCACATTTGGCATGGGACTTTCGGCCCTGTCCGCGTCCGTGTCGGCCAGACAAGCTTTGGGCATTGGCCACAATCGGGCCACAATCGAAAGCCGAGCAG",
"GGCAGCTGTCGGCAACTGTAAGCCATTTCTGGGACTTTGCTGTGAAAAGCTGGGCGATGGTTGTGGACCTGGACGAGCCACCCGTGCGATAGGTGAGATTCATTCTCGCCCTGACGGGTTGCGTCTGTCATCGGTCGATAAGGACTAACGGCCCTCAGGTGGGGACCAACGCCCCTGGGAGATAGCGGTCCCCGCCAGTAACGTACCGCTGAACCGACGGGATGTATCCGCCCCAGCGAAGGAGACGGCG",
"TCAGCACCATGACCGCCTGGCCACCAATCGCCCGTAACAAGCGGGACGTCCGCGACGACGCGTGCGCTAGCGCCGTGGCGGTGACAACGACCAGATATGGTCCGAGCACGCGGGCGAACCTCGTGTTCTGGCCTCGGCCAGTTGTGTAGAGCTCATCGCTGTCATCGAGCGATATCCGACCACTGATCCAAGTCGGGGGCTCTGGGGACCGAAGTCCCCGGGCTCGGAGCTATCGGACCTCACGATCACC"
]

# Call GreedyMotifSearch(Dna, k, t) and store the output in a variable called Motifs
Motifs = GreedyMotifSearch(Dna, 15, len(Dna))
# Print the Motifs variable and score
print(Motifs)
print(Score(Motifs))

#cool entropy function specifically for dna motifs
import math
def Entropy(motifs): #check coursera for actual mathematical explanation and theory
    entropy = 0
    for i in range(len(motifs)):
        for j in motifs[i]:
            if j == 0:
                entropy += 0
            else:
                entropy += j*(math.log(j,2))
    return entropy *-1


profile =[
[0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 0.9, 0.1, 0.1, 0.1, 0.3, 0.0],
[0.1, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.1, 0.2, 0.4, 0.6],
[0.0, 0.0, 1.0, 1.0, 0.9, 0.9, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0],
[0.7, 0.2, 0.0, 0.0, 0.1, 0.1, 0.0, 0.5, 0.8, 0.7, 0.3, 0.4]
]
print(Entropy(profile))

###################
#not part of course but here is RNA to protein code and then protein to RNA code:
def RNA_P(inp):
    ct = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    p = ""
    if len(inp)%3 == 0:
        for i in range(0, len(inp), 3):
            codon = inp[i:i + 3]
            p += ct[codon]
    return p

#def P_RNA(inp):
    #ct = {inverse codon table}
    #rna = ""
    #for i in range(0, len(inp)):

prof2 = [
[0.4, 0.3, 0.0, 0.1, 0.0, 0.9],
[0.2, 0.3, 0.0, 0.4, 0.0, 0.1],
[0.1, 0.3, 1.0, 0.1, 0.5, 0.0],
[0.3, 0.1, 0.0, 0.4, 0.5, 0.0]
]
#Pr("TCGGTA", prof2)

########################################################
# Continuation:

# Input:  A set of kmers Motifs (accounting Lapace's rule of sucession)
# Output: CountWithPseudocounts(Motifs)
def CountWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    count = {} # initializing the count dictionary
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(1)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count

# Input: A set of kmers Motifs
# Output: ProfileWithPseudocounts(Motifs)
def ProfileWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    cont=CountWithPseudocounts(Motifs)
    for i in range(k):
        su=0
        for symbol in "ACGT":
            su=su+cont[symbol][i]
        for symbol in "ACGT":
            cont[symbol][i] = cont[symbol][i]/su
    profile=cont
    return profile

# Here's my code that worked, please note in the ProfilewithMostProbablePattern, 
# you need to set the variable to -1 from 0, this is the bug that needs to be fixed in the code for it to work
def GreedyMotifSearchWithPseudocounts(Dna,k,t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for m in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][m:m+k])
        for j in range(1, t):
            P = ProfileWithPseudocounts(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

                    #################################################################
##################### RANDOM MOTIF (even more accurate! based partly on randomness) ################################
                    #################################################################

def Motifs(Profile, Dna):
    motifs = []
    t = len(Dna)
    k = 4 #can be customized I think
    for i in range(t):
        motif = ProfileMostProbablePattern(Dna[i], k, Profile)
        motifs.append(motif)
    return motifs
####################################################
def ProfileMostProbablePattern(Text, k, Profile):  #
    p_dict = {}                                    #
    for i in range(len(Text)- k +1):               #
        p = Pr(Text[i: i+k], Profile)              #
        p_dict[i] = p                              #
    m = max(p_dict.values())                       #
    keys = [k for k,v in p_dict.items() if v == m] #
    ind = keys[0]                                  #
    return Text[ind: ind +k]                       #
####################################################

#Generating random motifs for overall function later
#t isnt actually really needed)
def RandomMotifs(Dna, k, t):
    Motifs = []
    n = len(Dna[0])
    num = random.randint(0,n-k)
    for i in range (t):
        Motifs.append(Dna[i][num:num+k])
    return Motifs

#actual randomized motif search
def RandomizedMotifSearch(Dna, k, t):
    M = RandomMotifs(Dna, k, t)
    BestMotifs = M
    while True:
        Profile = ProfileWithPseudocounts(M)
        M = Motifs(Profile, Dna)
        if Score(M) < Score(BestMotifs):
            BestMotifs = M
        else:
            return BestMotifs

                    ###################################################################
##################### Gibbs sampling method for finding hidden motifs: the final boss ################################
                    ###################################################################
                
#first prerequisite, to transform a dictionary of pseudocount probabilities evenly so they sum to 1 always.
# Input: A dictionary Probabilities, where keys are k-mers and values are the probabilities of these k-mers (which do not necessarily sum up to 1)
# Output: A normalized dictionary where the probability of each k-mer was divided by the sum of all k-mers' probabilities
def Normalize(P):
    D = {}
    for k,v in P.items():
        D[k] = P[k]/sum(P.values())
    return D

#example with 5, but usually U shouldnt be there d = {'A': 0.45, 'C': 0.63, 'G': 0.09, 'T':0.27, 'U':0.36}

# Input:  A dictionary Probabilities (That have been normalized by previous function!!!) whose keys are k-mers and whose values are the probabilities of these kmers
# Output: A randomly chosen k-mer with respect to the values in Probabilities
def WeightedDie(Probabilities):
    import random
    randomn = random.uniform(0,1)
    for p in Probabilities:
        randomn-= Probabilities[p]
        if randomn <= 0:
            return p

# Input:  A string Text, a profile matrix Profile, and an integer k
# Output: ProfileGeneratedString(Text, profile, k)
def ProfileGeneratedString(Text, profile, k):
    n = len(Text)
    probabilities = {}
    for i in range(0,n-k+1):
        probabilities[Text[i:i+k]] = Pr(Text[i:i+k], profile)
    probabilities = Normalize(probabilities)
    return WeightedDie(probabilities)

#Be aware, the the gibbs algorithm has the common problem of running into a potential local minimum instead of global minimum
# it should therefore be run many times.
#alternatively, a more advanced version could be implemented with packages allowing calculus, 
#and extrapolation of functions through data points with basic ML.

#Here is the function from the course:
def GibbsSampler(Dna, k, t, N):
    BestMotifs = [] # output variable
    Motifs = RandomMotifs(Dna, k, t)
    #print(Motifs)
    BestMotifs = Motifs
    for j in range(N):
        i = random.randint(0,t-1)
        new_Motif = []
        for k1 in range(t):
            if k1!=i:
                new_Motif.append(Motifs[k1])
        profile = ProfileWithPseudocounts(new_Motif)
        motif_i = ProfileGeneratedString(Dna[i], profile, k)
        Motifs[i] = motif_i
        if Score(Motifs)<Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs
