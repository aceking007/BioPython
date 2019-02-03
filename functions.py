import os
import re
from math import log

## Defining the protein tuple
protein = ('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')
##protein = set(protein)

## Defining the DNA tuple
dna = ('A', 'T', 'G', 'C')
##dna = set(dna)

## Defining the RNA tuple
rna = ('A', 'U', 'G', 'C')
##rna = set(rna)

## Saving the genetic code to a dictionary
## Requisites: gencode file
dnapro = {}
fd = open('gencode', 'r')
fd.seek(0)
gencode = fd.readlines()
fd.close()
for i in range(len(gencode)):
    gencode[i] = gencode[i].strip()
    c = gencode[i].split()
    dnapro[c[0]] = c[1]

## Function to test if the sequence is DNA or RNA or Protein
## Input: Sequence
## Output: dna, rna, protein or null
def checkSeq(i):
    chk = set(i)
    d = set(dna)
    r = set(rna)
    p = set(protein)
    if chk.intersection(d) == chk:
        return 'dna'
    elif chk.intersection(r) == chk:
        return 'rna'
    elif chk.intersection(p) == chk:
        return 'protein'
    else:
        return 'null'

## Function to return length of sequence
## Input: Sequence
## Output: Length
def lenSeq(i):
    return len(i)

## Function to return occurence of each nucleotide or amino acid
## Input: Sequence
## Output: Dictionary with all units and occurences
def occur(i):
    dict = {}
    if checkSeq(i) == 'dna':
        dict = {dna[x]:0 for x in range(len(dna))}
        for l in i:
            if l in dict:
                dict[l] += 1
    elif checkSeq(i) == 'rna':
        dict = {rna[x]:0 for x in range(len(rna))}
        for l in i:
            if l in dict:
                dict[l] += 1
    elif checkSeq(i) == 'protein':
        dict = {protein[x]:0 for x in range(len(protein))}
        for l in i:
            if l in dict:
                dict[l] += 1
    else:
        return {}
    return dict

## Function to calculate AT/GC value
## Input: DNA Sequence
## Output: 0 or AT/ATGC
def atcontent(i):
    if checkSeq(i) == 'dna':
        a = occur(i)
        b = float(a['A'] + a['T'])/float(lenSeq(i))
        return b
    else:
        return 0

## Function for DNA complement
## Input: DNA Sequence
## Output: Complementary DNA Sequence
def cmplmnt(a):
    l = ''
    comd = {'A':'T','T':'A','G':'C','C':'G'}
    if checkSeq(a) == 'dna':
        for i in a:
            l += comd[i]
        return l
    else:
        print 'It is not a DNA sequence'

## Function to find occurence of a given motif
## Input: DNA or RNA Sequence
## Output: Number of matches & Span of all matches
def matchSeq(i):
    if checkSeq(i) == 'dna' or checkSeq(i) == 'rna':
        m = raw_input('Enter the motif you want to find:')
        a = re.finditer(m,i,re.I)
        count = 0
        for i in a:
            print i.span()
            count += 1
        print 'There were ' + str(count) + ' matches for the given motif in the sequence'
    else:
        return 0

## Function to represent part of sequence in lowercase
## Input: Sequence and part for processing
## Output: Modified sequence
def lowSeq(a):
    if checkSeq(a) == 'dna' or checkSeq(a) == 'rna' or checkSeq(a) == 'protein':
        m = raw_input('Enter the sequence to represent in lowercase:')
        b = re.sub(m,m.lower(),a,flags=re.IGNORECASE)
        print b
    else:
        print 'This is an incorrect sequence'

## Function to calculate Shannon entropy
## Input: Sequence
## Output: Shannon Entropy
def shannonEntropy(a):
    dict = occur(a)
    h = 0
    if bool(dict):
        b = lenSeq(a)
        for i in dict:
            val = dict[i]
            dict[i] = float(val)/float(b)
        for i in dict:
            if dict[i] != 0:
                h += dict[i]*(log(dict[i],2))
        h = h*(-1)
        return h
    else:
        return 'This is not a valid sequence'

## Function to extract specific information from FASTA dictionary
## Input type: dict({Description:Sequence})
## Output: list[(Seq_name, Accession_no, Organism, Seq)]
def fasta(a):
    l = []
    names = []
    acs = []
    org = []
    sq = []
    for i in a:
        l.append(i)
        sq.append(a[i])
        b = l[0].split('|')
        acs.append(b[1])
        m = re.compile(r'\w+')
        name = m.search(b[2]).group()
        names.append(name)
        n = re.compile(r'.*(?<=OS=)(.+)(?=\s+OX)')
        org.append(n.search(b[2]).group(1))
        l = []
    mapped = zip(names, acs, org, sq)
    return mapped

## Function to divide string into small chunks
## Input: s == string, w == size of each chunk
## Output: List
def wrap(s, w):
    return [s[i:i+w] for i in range(0, len(s), w)]


## Function to generate expected Protein from DNA
## Input: dict({Description, Sequence})
## Prints Sequence description, Sequence, Most probable protein
def dnatopro(d):
    l = ''
    for i in d:
        if checkSeq(d[i]) == 'dna':
            print ''
            print 'Sequence description: ' + i
            print 'DNA Sequence: ' + d[i]
            m = re.search(r'ATG',d[i]).start()
            part = wrap(d[i][m:],3)
            for i in part:
                if i == 'ATG':
                    count = 1
                if i in dnapro and count == 1:
                    if dnapro[i] == '_':
                        break
                    else:
                        l += dnapro[i]
            print ''
            print 'Generated protein: ' + l
            print ''
            l = ''
        else:
            print i
            print 'This is not a valid DNA sequence'
            print ''

## Function to read Blosum Matrix
## Output: 2-D Dictionary containing all the values
def blosum():
    file = open('blosum','r')
    file.seek(0)
    l = []
    for i in file:
        if i[0] != '#':
            l.append(i)
    file.close()
    for i in range(len(l)):
        l[i] = l[i].strip()
    a = l[0].split()
    mat = {}
    for i in range(1,len(l)):
        b = l[i].split()
        for j in range(1,len(b)):
            b[j] = int(b[j])
        mapped = dict(zip(a,b[1:]))
        mat[a[i-1]] = mapped
    return mat
