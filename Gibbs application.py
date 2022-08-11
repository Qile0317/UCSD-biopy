# import the random package
import random
# Copy your GibbsSampler function (along with all required subroutines) below this line
def GibbsSampler(Dna, k, t, N):
    BestMotifs = [] # output variable
    # your code here
    Motifs = RandomMotifs(Dna, k, t)
    BestMotifs = Motifs
    for j in range(1,N):
        i = random.randint(0,t-1)
        ReducedMotifs = []
        for j in range(0,t):
            if j != i:
                ReducedMotifs.append(Motifs[j])
        Profile = ProfileWithPseudocounts(ReducedMotifs)
        Motif_i = ProfileGeneratedString(Dna[i], Profile, k)
        Motifs[i] = Motif_i
        if Score(Motifs) < Score(BestMotifs):
                BestMotifs=Motifs
    return BestMotifs

# place all subroutines needed for GibbsSampler below this line
# Input:  A list of strings Dna, and integers k and t
# Output: RandomMotifs(Dna, k, t)
# HINT:   You might not actually need to use t since t = len(Dna), but you may find it convenient

def RandomMotifs(Dna, k, t):
    # place your code here.
    s = len(Dna[0])
    rm = []
    for i in range(0,t):
        init_index = random.randint(1,s-k)
        rm.append(Dna[i][init_index:init_index+k])    
    return rm

# Input:  A set of kmers Motifs
# Output: ProfileWithPseudocounts(Motifs)

def ProfileWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {} # output variable
    # your code here
    c = CountWithPseudocounts(Motifs)
    for n in 'ACGT':
        p = []
        for i in range(0,k):
            p.append(c[n][i]/(t+4))
        profile[n] = p
    return profile

# Input:  A set of kmers Motifs
# Output: CountWithPseudocounts(Motifs)

def CountWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    # insert your code here
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

#tests in which of the intervals defined by list ar the number r lies

def testinterval(ar,r):
    ar.sort()
    if r<= ar[0]:
      return ar[0]
    for i in range(1,len(ar)-1):
      if ar[i-1]<r<=ar[i]:
        return ar[i]
    if ar[len(ar)-2]< r:
      return ar[len(ar)-1]

# Input:  A dictionary Probabilities whose keys are k-mers and whose values are the probabilities of these kmers
# Output: A randomly chosen k-mer with respect to the values in Probabilities

def WeightedDie(Probabilities):
    # your code here
    sumprob = {}
    s = 0
    for p in Probabilities:
        s += Probabilities[p]
        sumprob[p] = s
    revprob = {}
    for q in sumprob:
      revprob[sumprob[q]] = q
    w = list(sumprob.values())
    r = random.uniform(0,1)
    kmer = revprob[testinterval(w,r)]
    return kmer

# Input:  A string Text, a profile matrix Profile, and an integer k
# Output: ProfileGeneratedString(Text, profile, k)

def ProfileGeneratedString(Text, profile, k):
    # your code here
    n = len(Text)
    probabilities = {} 
    for i in range(0,n-k+1):
        probabilities[Text[i:i+k]] = Pr(Text[i:i+k], profile)
    probabilities = Normalize(probabilities)
    return WeightedDie(probabilities)

# Input:  String Text and profile matrix Profile
# Output: Pr(Text, Profile)

def Pr(Text, Profile):
    # insert your code here
    p = 1
    for i in range(0,len(Text)):
        p *= Profile[Text[i]][i]
    return p

# Input: A dictionary Probabilities, where keys are k-mers and values are the probabilities of these k-mers (which do not necessarily sum up to 1)
# Output: A normalized dictionary where the probability of each k-mer was divided by the sum of all k-mers' probabilities

def Normalize(Probabilities):
    # your code here
    result = {}
    sum = 0
    for m in Probabilities:
        sum += Probabilities[m]
    for n in Probabilities:
        result[n]= Probabilities[n]/sum
    return result  

# Input:  A set of k-mers Motifs
# Output: The score of these k-mers.

def Score(Motifs):
    # Insert code here
    k = len(Motifs[0])
    t = len(Motifs)
    cs = ConsensusWithPseudocounts(Motifs)
    score = 0
    for j in range(0,k):
        for i in range(0,t):
            if Motifs[i][j] != cs[j]:
                score += 1
    return score

# Input:  A set of kmers Motifs
# Output: A consensus string of Motifs.

def ConsensusWithPseudocounts(Motifs):
    # insert your code here
    k = len(Motifs[0])
    count = CountWithPseudocounts(Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus    
k = 8
t = 5
N = 100
Dna = ["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA","GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG","TAGTACCGAGACCGAAAGAAGTATACAGGCGT","TAGATCAAGTTTCAGGTGCACGTCGGTGAACC","AATCCACCAGCTCCACGTGCAATGTTGGCCTA"]

# Copy the ten strings occurring in the hyperlinked DosR dataset below.
Dna = ["GCGCCCCGCCCGGACAGCCATGCGCTAACCCTGGCTTCGATGGCGCCGGCTCAGTTAGGGCCGGAAGTCCCCAATGTGGCAGACCTTTCGCCCCTGGCGGACGAATGACCCCAGTGGCCGGGACTTCAGGCCCTATCGGAGGGCTCCGGCGCGGTGGTCGGATTTGTCTGTGGAGGTTACACCCCAATCGCAAGGATGCATTATGACCAGCGAGCTGAGCCTGGTCGCCACTGGAAAGGGGAGCAACATC", 
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

# set t equal to the number of strings in Dna, k equal to 15, and N equal to 100

t = 10
k = 15
N = 100
# Call GibbsSampler(Dna, k, t, N) 20 times and store the best output in a variable called BestMotifs

BestMotifs = GibbsSampler(Dna, k, t, N)

# Print the BestMotifs variable

print(BestMotifs)

# Print Score(BestMotifs)

print(Score(BestMotifs))