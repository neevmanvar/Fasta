f = open('dna2.fasta','r')  #open a fasta file

seq = f.read()              # read sequence
record = seq.count('>')     # number of total sequences
print("record is",record)

# Build a dictionary with ID
seqs={}
f.seek(0)
for line in f:
    line = line.rstrip()    # right side split
    if line[0]=='':         #skip any empty line
        continue
    elif line[0]=='>':
        words = line.split()
        name = words[0][1:]
        seqs[name]=''
    else:
        seqs[name] = seqs[name]+line


list = []
new_list = []
for name,a in seqs.items():
    new_list.append(len(a))
    list.append(name)

print("Final length of all sequences:",new_list)


max_seq = max(new_list)
min_seq = min(new_list)

for index,Max in enumerate(new_list):
    if Max == max_seq:
        break
maximum = ([list[index],max_seq])

for index,Min in enumerate(new_list):
    if Min == min_seq:
        break
minimum = ([list[index],min_seq])
print("Maximum sequence and it's ID is",maximum)
print("Minimum sequence and it's ID is",minimum)

## Finding ORF in all three frames

def has_stop_codon(dna,frame):
    stop_codon_found = False
    stop_codons = ['TGA','TUG','TAA']

    for i in range(frame,len(dna),3):
        codon = dna[i:i+3]
        if codon in stop_codons:
            stop_codon_found = True
            break
    return stop_codon_found

def has_start_codon(dna,frame):
    start_codon_found = False
    start_codon = ['ATG']

    for i in range(frame,len(dna),3):
        codon = dna[i:i+3]
        if codon in start_codon:
            start_codon_found = True
            break
    return start_codon_found



def ORF_find(dna,frame):

    stop_codons = ['TGA','TUG','TAA']
    start_codon = ['ATG']
    ORF_list = []
    list_ORF_length = []
    list_ORF_index = []
    for i in range(frame,len(dna),3):
        codon = dna[i:i+3]                  # slicing dna
        if codon in start_codon:            #checking start codon
            for j in range(i+3,len(dna),3):    # start from nect sequence
                codon2 = dna[j:j+3]             # second sequencing
                if codon2 in stop_codons:       #checking stop codons
                    ORF = dna[i:j+3]            # total ORF
                    ORF_length = len(ORF)       # ORF length
                    ORF_list.append(ORF)        # append on list
                    list_ORF_length.append(ORF_length)  #append ORF lengths on list
                    list_ORF_index.append(dna.find(ORF))   #append ORF  posion on list
                    break
    return ORF_list,list_ORF_length,list_ORF_index


# code for getiing maximum and minimum ORF in given dna sequence
ORF_max = {}
ORFs = {}
ORF_min = {}
ORF_index_max = {}
ORF_index_min = {}
ORF_numbers = {}

for name,dna in seqs.items():
    frame =2
    f = str(frame)
    if (has_start_codon(dna,frame) and has_stop_codon(dna,frame)):  # checking start and stop codons in dna sequence
        ORF_list,list_ORF_length,list_ORF_index = ORF_find(dna,frame)              # find ORF
        if (list_ORF_length != [] and list_ORF_index != []):                                   # check if ORF length is not zero
            ORF_max[name] = max(list_ORF_length)  # use max to get maximum length
            ORF_min[name] = min(list_ORF_length)  # use min to get minimum length

            max_index = list_ORF_length.index(max(list_ORF_length))
            ORF_index_max[name] = list_ORF_index[max_index]  #find exact positions of max ORFs

            min_index = list_ORF_length.index(min(list_ORF_length))
            ORF_index_min[name] = list_ORF_index[min_index]  #find exact position of min ORFs

        ORF_numbers[name] = len(ORF_list)         # using len fun to get numbers of ORF in given sequence


# Find max ORF with ID and position
init = 0
index = 0
for name,length in ORF_max.items():
    if init < length:
        init = length
        new = name
        max_position = ORF_index_max[new] + 1

    index += 1

print('\n In frame ' + f,'Longest ORFs ID: \t',new,'\n','Longest ORFs length: \t',init,'\n','Positon of Longest ORF: \t',max_position)


# FInd min ORF with ID and position

init = float('inf')
index = 0
for name,length in ORF_min.items():
    if init > length:
        init = length
        new = name
        max_position = ORF_index_max[new] + 1

    index += 1

print('\n In frame ' + f,'shortest ORFs ID: \t',new,'\n','shortest ORFs length: \t',init,'\n','Positon of shortest ORF: \t',max_position)

#  What is the length of the longest forward ORF that appears in the sequence with the identifier gi|142022655|gb|EQ086233.1|16

print('length of the longest forward ORF that appears in the sequence with the identifier is: \t',ORF_max['gi|142022655|gb|EQ086233.1|16'])


# Find repeats in the dna sequence with length of finit number


def find_repeat(dna,L):
    Repeats = {}

    for i in range(0,len(dna)):
        if dna[i:i+L] not in Repeats:
            Repeats[dna[i:i+L]] = 1
        else:
            Repeats[dna[i:i+L]] += 1
    return Repeats

# create dictionary of most freq word
Repeat_words = {}
World_length = []
D = {}

L = 6       # dna sequence with finit length
Total_L = str(L)

for name,dna in seqs.items():
    max_value = 0
    Repeats = find_repeat(dna,L)
    for key,value in Repeats.items():
        if key not in D:        #store same sequences in another dictionary
            D[key] = value
        else:
            D[key] += value     #add new value to the old value

for key,value in D.items():         #unpack dictionary
    if max_value < value and len(key) == L:     #max value condition
        max_value = value
        new_key = key


print('the most frequent word with length of '+ Total_L+' is": \t',new_key,' ',max_value)


# In single dna sequence highest number of length occure

for name,dna in seqs.items():
    max_value = 0
    Repeats = find_repeat(dna,L)
    for key,value in Repeats.items():
        if max_value < value and len(key) == L:
            max_value = value
            new_key = key
    Repeat_words[name] = [new_key,max_value]

print('Maximum sequence of length '+Total_L+' occure in each dna sequence:\n',Repeat_words)
