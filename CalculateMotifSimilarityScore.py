"""

This script scans the human TF motifs (provided as position weight matrix PWM) for a given sequence
The output is the similarity score for each sliding window with step set as 1bp

@author: Sunny Sun. email: sunny.sun@nyu.edu

To Use: 
Create two folders in the same directory: Promoter_fasta, pwm
Promoter_fasta contains all the input sequences
pwm contains all the motif matrix

Running the script:
python CalculateMotifSimilarityScore.py -i [name of the input sequence (including full path)]

"""

################################################################################
################################ Modules #######################################
################################################################################

import re
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import pandas as pd
from os import listdir, makedirs
from os.path import isfile, join, exists
import optparse

################################################################################
################################ Functions #####################################
################################################################################


# Extracts the sequence in a FASTA file.
def load_FASTA (filename):
    """
    This function takes a filename (str) as input and returns
    the sequence (str) found in that file. 
    """
    # open file according to the filename given.
    with open (filename) as f:
       # reads what's in the file as a string, and then splits
       # the string by '>' into a list.
       lines = f.read().split('>')[1:]
       # close the file.
       f.close()
       # print what's after the first '\n'
       return ''.join(lines[0].split('\n')[1:])


# Function to read in the pwm
def load_pwm (motif_pwm_filename):
    # read in the pwm file, each position is the probability
    pwm = pd.read_csv(motif_pwm_filename,sep='\t',header=None,engine='python')
    pwm_len=len(pwm)
    # insert the 5th column for N
    pwm.insert(4,4,0)
    # rename the columns
    pwm.columns = ['A', 'C', 'G', 'T', 'N']
    # create a dictionary for 5 letters
    pwm_order={'A':0,'C':1,'G':2,'T':3,'N':4}
    # return the matrix, length and dictionary
    return pwm,pwm_len,pwm_order


# calcualte the similarity score for a given sequence in each bp
# each sliding window, minus the first nucleotide, plus the last nucleotide
# width of the sliding window = length of motif
def similarity_score (sequence, pwm_len, pwm, pwm_order, pwm_dic, pwm_max_sum):
    # add this line if input is not in the Seq format
    inp=Seq(sequence, generic_dna)
    # change the letters to the upper
    inp=inp.upper()
    
    # creates a dataframe to store output scores
    out_len=len(inp)-pwm_len+1
    pos=list(range(out_len))
    df=pd.DataFrame()
    df['pos']=pos
    score=[]
    score_rc=[]
    
    # calculate each score in a pwm_len window, step is 1bp
    for k in range(out_len):
        seq=inp[k:(k+pwm_len)]
        # find each bp probablity using dictionary
        keys=[str(m)+n for m,n in zip(range(pwm_len),[i for i in seq])]
        fe=[pwm_dic.get(key) for key in keys]
        # sum up the probability of each position
        addup=sum(fe)
        score.append(addup)      
    df['score']=[i/pwm_max_sum for i in score]
    
    ## calculate scores of the reverse complement
    inp_rc=inp.reverse_complement()
    # calculate each score in a pwm_len window, step is 1bp
    for k in range(out_len):
        seq=inp_rc[k:(k+pwm_len)]
        # find each bp probablity using dictionary
        keys=[str(m)+n for m,n in zip(range(pwm_len),[i for i in seq])]
        fe=[pwm_dic.get(key) for key in keys]
        # sum up the probability of each position
        addup=sum(fe)
        score_rc.append(addup)
    
    # for each position, compare the scores on both strands, take the higher one
    score_rc.reverse()
    df['score_rc']=[i/pwm_max_sum for i in score_rc]
    # return the max ratio for each position
    df['similarityScore']=[0]*out_len
    for index, row in df.iterrows():
        df.loc[index,'similarityScore']=max(row['score'],row['score_rc'])
    df=df.drop(['score'],axis=1)
    df=df.drop(['score_rc'],axis=1)
    df['pos']=[i+1 for i in pos]
    # output a dataframe with two columns: position, score
    return df


# loop through all the motifs for a input sequence
def motif_matches(input_seq_filename):
    """
    This function takes one input sequence, calculates all the motif match scores at each position and 
    for all the motifs.
    input_seq_filename needs to include the full path
    """
    # read in the input sequence
    seq = load_FASTA(input_seq_filename)

    # creates all the folders to store results
    #basepath='/Users/sunnysun/Desktop/wenfeng/'
    # finds the basename, make sure it ends with '/'
    basepath=input_seq_filename.split('Promoter_fasta')[0]
    if not basepath.endswith('/'):
        basepath=basepath+'/'
    # make a new folder motifMatches to store all the scores for each sequence
    if not exists(basepath+'motifMatches'):
        makedirs(basepath+'motifMatches')
    # make a new folder MaxScoreAndMatches to store the output summary table
    if not exists(basepath+'MaxScoreAndMatches'):
        makedirs(basepath+'MaxScoreAndMatches')
    # creates a folder for each input sequence
    newpath = basepath+'motifMatches/' + input_seq_filename.split('/')[-1].split('.')[0]
    if not exists(newpath):
        makedirs(newpath)
    # creates a folder to store all the motif scores
    newpath_allscores=newpath+'/Allscores'
    if not exists(newpath_allscores):
        makedirs(newpath_allscores)
    # creates a folder to store all the matches sits
    newpath_matches=newpath+'/MatchedSites'
    if not exists(newpath_matches):
        makedirs(newpath_matches)
    
    # list all the motif names in the pwm
    motifpath=basepath+'pwm/'
    motiffiles = [f for f in listdir(motifpath) if isfile(join(motifpath, f))]
    
    # create a summary table with three columns
    allmotifs=pd.DataFrame()
    allmotifs['TF']=motiffiles
    MaxScore=[]
    NumberOfMatches90=[]
    #NumberOfMatches85=[]
    
    # Loop through all the motifs, for each one, calculate the similarity score at each position
    for j in motiffiles:
        # calculate the maximum of summed weights
        pwm=load_pwm(motifpath+j)[0]
        pwm_len=load_pwm(motifpath+j)[1]
        pwm_order=load_pwm(motifpath+j)[2]
        
        # create a dictionary for pwm
        pwm_dic={}
        for i in range(pwm_len):
            pwm_dic.update({(str(i)+'A'): pwm.loc[i,'A']})
            pwm_dic.update({str(i)+'C': pwm.loc[i,'C']})
            pwm_dic.update({str(i)+'G': pwm.loc[i,'G']})
            pwm_dic.update({str(i)+'T': pwm.loc[i,'T']})
            pwm_dic.update({str(i)+'N': pwm.loc[i,'N']})
        pwm_max_sum=0
        
        # calculate the maximum score for each motif
        for m in range(pwm_len):
            pwm_max_sum+=(max(pwm.loc[m,:]))
        score_matrix=similarity_score(seq, pwm_len, pwm, pwm_order, pwm_dic, pwm_max_sum)
        
        # set the cutoff with 0.9 for the similarity score
        score_matrix_90=score_matrix.loc[score_matrix['similarityScore']>=0.9,]
        #score_matrix_85=score_matrix.loc[score_matrix['similarityScore']>=0.85,]
        
        # calculate the maximum score
        max_score=max(score_matrix['similarityScore'])
        MaxScore.append(max_score)
        
        # calcualte the number of matched sites
        numberofSites90=score_matrix_90.shape[0]
        NumberOfMatches90.append(numberofSites90)
        #numberofSites85=score_matrix_85.shape[0]
        #NumberOfMatches85.append(numberofSites85)
        
        # write out the similarity score at each position
        score_matrix.to_csv(newpath_allscores+'/'+input_seq_filename.split('/')[-1].split('.')[0]+'-'+j+'-MotifScores.txt',sep='\t',header=True,index=False)
        # write out the matched position
        score_matrix_90.to_csv(newpath_matches+'/'+input_seq_filename.split('/')[-1].split('.')[0]+'-'+j+'-MotifScores90%.txt',sep='\t',header=True,index=False)
        #score_matrix_85.to_csv(newpath+'/'+input_seq_filename.split('/')[-1].split('.')[0]+'-'+j+'-MotifScores85%.txt',sep='\t',header=True,index=False)
    
    # outputs a table with three columns: TF_name, MaxScore, NumberOfMatches
    allmotifs['MaxScore']=MaxScore
    allmotifs['NumberOfMatches90%']=NumberOfMatches90
    #allmotifs['NumberOfMatches85%']=NumberOfMatches85
    allmotifs.to_csv(basepath+'MaxScoreAndMatches/'+input_seq_filename.split('/')[-1].split('.')[0]+'-MaxScoreAndMatches.txt',sep='\t',header=True,index=False)
   

################################################################################
################################# MAIN #########################################
################################################################################
# main

# parse object for managing input options.      
parser = optparse.OptionParser()

# essential data, defines a commanline option "-i"
parser.add_option('-i', dest = 'input_seq_filename', default = '', help = 'This input\
 is the name of the input sequence filename (including full path)')

# loads the inputs
(options, args) = parser.parse_args()

# reads the inputs from command lines
input_seq_filename = options.input_seq_filename

# runs the function
motif_matches(input_seq_filename)




