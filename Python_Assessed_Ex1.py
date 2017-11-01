# set up dictionaries and lists to store information
threshold = 30      # Q score threshold
seq_pair = {}       # match sample number to both read pairs
sequence = {}       # match sample number to the sequence
quality = {}        # match sample number to the ascii values
Qscore = {}         # match sample number to calculated averaged Q score
seq_Qscore30_more = []      # list of sample number with Q score of both read pairs more than 30
seq_Qscore30_less = []      # list of sample number with Q score of at least one of the read pairs less than 30

# extraction of information from the original file
with open('seq_sample.fastq') as in_file:
    raw_file = in_file.readlines()      # open file and read lines
    for count in range(0,len(raw_file),8):      # every 8 lines represents all of the information of one sample
        seq_pair[raw_file[count].rstrip()[:-2]] = [raw_file[count].rstrip(), raw_file[count+4].rstrip()]
        # match the sample number with both read pairs name
        sequence[raw_file[count].rstrip()[:-2]] = [raw_file[count + 1].rstrip(), raw_file[count + 5].rstrip()]
        # match the sample number with both sequence of read pairs
        quality[raw_file[count].rstrip()[:-2]] = [raw_file[count+3].rstrip(), raw_file[count+7].rstrip()]
        # match the sample number with both ascii value of read pairs
in_file.close()

# calculation of the Q score
for key in sequence:
    # select of the ascii quality score of the pair reads separately
    qual1 = quality[key][0]
    qual2 = quality[key][1]
    sum1 = sum2 = 0     # initiate summation of the Q scores

    for count1 in range(0, len(qual1)):
        sum1 += ord(qual1[count1]) - 64     # converting the ascii score to Q score and adding up all the Q scores
    for count2 in range(0, len(qual2)):
        sum2 += ord(qual2[count2]) - 64     # converting the ascii score to Q score and adding up all the Q scores

    Qscore[key] = [sum1/len(qual1), sum2/len(qual2)]
    # calculate average Q score, store them to the dictionary in a list

    if sum1/len(qual1) >= threshold and sum2/len(qual2) >= threshold:
        # if both of the pair reads Q score above or equal to the threshold then store the sample name to the list
        seq_Qscore30_more.append(key)
    else:
        # otherwise (at least one of the pair read Q score less than 30) store to seq_Qscore30_less
        seq_Qscore30_less.append(key)

# output the information to two separate files using out_file
out_file = open("seq_sample_Qscore30+.fastq", "w")
for key in seq_Qscore30_more:       # elements in seq_Qscore30_more are the sample names with Q score>=30
    for i in [0,1]:     # for loop to get both of the pair reads
        # following the sample format as that in seq_sample.fastq
        out_file.write(seq_pair[key][i]+'\n')       # start with the pair read names
        out_file.write(sequence[key][i] + '\n')     # sequence of the pair read
        out_file.write('+'+'\n')                    # '+' sign
        out_file.write(quality[key][i]+'\n')        # quality in ascii values
out_file.close()        # close the file

# repeat for the ones that at least one of the pair read Q score less than 30
out_file = open("seq_sample_Qscore30-.fastq", "w")
for key in seq_Qscore30_less:
    for i in [0,1]:
        out_file.write(seq_pair[key][i]+'\n')
        out_file.write(sequence[key][i] + '\n')
        out_file.write('+'+'\n')
        out_file.write(quality[key][i]+'\n')
out_file.close()