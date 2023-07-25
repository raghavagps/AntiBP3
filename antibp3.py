######################################################################################
# AntiBP3 is developed for predicting , Desigining and scanning of Antibacterial peptides  #
# It is developed by Prof G. P. S. Raghava's group.                #
# Please cite: https://webs.iiitd.edu.in/raghava/antibp3/                        #
######################################################################################
import argparse
import warnings
#import pkg_resources
import os
import sys
import numpy as np
import pandas as pd
import itertools
import pickle
import re
import glob
import time
import uuid
import zipfile
warnings.filterwarnings('ignore')

# Defining the function for generating all possible mutants
def mutants(file1,file2):
    std = list("ACDEFGHIKLMNPQRSTVWY")
    cc = []
    dd = []
    ee = []
    df2 = file2
    df2.columns = ['Name']
    df1 = file1
    df1.columns = ['Seq']
    for k in range(len(df1)):
        cc.append(df1['Seq'][k])
        dd.append('Original_'+'Seq'+str(k+1))
        ee.append(df2['Name'][k])
        for i in range(0,len(df1['Seq'][k])):
            for j in std:
                if df1['Seq'][k][i]!=j:
                    dd.append('Mutant_'+df1['Seq'][k][i]+str(i+1)+j+'_Seq'+str(k+1))
                    cc.append(df1['Seq'][k][:i] + j + df1['Seq'][k][i + 1:])
                    ee.append(df2['Name'][k])
    xx = pd.concat([pd.DataFrame(ee),pd.DataFrame(dd),pd.DataFrame(cc)],axis=1)
    xx.columns = ['Seq ID','Mutant ID','Seq']
    return xx

# defining the function to check the seqeunce
def readseq(file):
    with open(file) as f:
        records = f.read()
    records = records.split('>')[1:]
    seqid = []
    seq = []
    for fasta in records:
        array = fasta.split('\n')
        name, sequence = array[0].split()[0], re.sub('[^ACDEFGHIKLMNPQRSTVWY-]', '', ''.join(array[1:]).upper())
        seqid.append('>'+name)
        seq.append(sequence)
    if len(seqid) == 0:
        f=open(file,"r")
        data1 = f.readlines()
        for each in data1:
            seq.append(each.replace('\n',''))
        for i in range (1,len(seq)+1):
            seqid.append(">Seq_"+str(i))
    df1 = pd.DataFrame(seqid)
    df2 = pd.DataFrame(seq)
    return df1,df2

# defining the function to check the length of seqeunces
def lenchk(file1):
    cc = []
    df1 = file1
    df1.columns = ['seq']
    for i in range(len(df1)):
        if len(df1['seq'][i])>50:
            cc.append(df1['seq'][i][0:50])
        else:
            cc.append(df1['seq'][i])
    df2 = pd.DataFrame(cc)
    df2.columns = ['Seq']
    return df2

def nct(file,n):
    df1 = file
    df2 = pd.DataFrame(df1.iloc[:,0].str.upper())
    df3 = []
    for i in range(0,len(df2)):
        df3.append(df2.iloc[i,0][:n]+df2.iloc[i,0][-n:][::-1])
        df4 = pd.DataFrame(df3)
    for i in range(0,len(df4)):
        ss = len(df4.iloc[:,0][i])
        if ss/2 < n:
            print('\nSequence number',i+1,'has length of',int(ss/2),'which is less than the provided value of N- and C-terminal, that is',n,'. Kindly check the sequences.')
            sys.exit()
    return df4

def aab(file,out):
    std = list("ACDEFGHIKLMNPQRSTVWY")
    df = pd.read_csv(file, header = None)
    uu = []
    for ss in df[0]:
         uu.append(len(ss))
    zz = df.iloc[:,0]
    f = open(out, mode='w')
    sys.stdout = f
    A=('1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0')
    C=('0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0')
    D=('0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0')
    E=('0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0')
    F=('0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0')
    G=('0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0')
    H=('0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0')
    I=('0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0')
    K=('0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0')
    L=('0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0')
    M=('0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0')
    N=('0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0')
    P=('0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0')
    Q=('0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0')
    R=('0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0')
    S=('0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0')
    T=('0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0')
    V=('0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0')
    W=('0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0')
    Y=('0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1')
    for mm in range (1,max(uu)+1):
        for ee in std:
            print(ee+str(mm),end=',')
    print("")
    for i in range(0,len(zz)):
        for j in zz[i]:
            if j == "A":
                print(''.join(A), end = ',')
            if j == "C":
                print(''.join(C), end = ',')
            if j == "D":
                print(''.join(D), end = ',')
            if j == "E":
                print(''.join(E), end = ',')
            if j == "F":
                print(''.join(F), end = ',')
            if j == "G":
                print(''.join(G), end = ',')
            if j == "H":
                print(''.join(H), end = ',')
            if j == "I":
                print(''.join(I), end = ',')
            if j == "K":
                print(''.join(K), end = ',')
            if j == "L":
                print(''.join(L), end = ',')
            if j == "M":
                print(''.join(M), end = ',')
            if j == "N":
                print(''.join(N), end = ',')
            if j == "P":
                print(''.join(P), end = ',')
            if j == "Q":
                print(''.join(Q), end = ',')
            if j == "R":
                print(''.join(R), end = ',')
            if j == "S":
                print(''.join(S), end = ',')
            if j == "T":
                print(''.join(T), end = ',')
            if j == "V":
                print(''.join(V), end = ',')
            if j == "W":
                print(''.join(W), end = ',')
            if j == "Y":
                print(''.join(Y), end = ',')
        print("")
    f.close()
    sys.stdout = sys.__stdout__

# defining the function to generate the features out of sequences
def feature_gen(file, ncter, wd):
    file1 = nct(file,ncter)
    file1.to_csv(wd + '/sam_input.csv', index=None, header=False)
    aab(wd + '/sam_input.csv', wd + '/tempfile_out')
    df = pd.read_csv(wd + '/tempfile_out')
    df.columns = 'NC'+df.columns
    df2 = df.iloc[:,:-1].fillna('NA')
    return df2


# Function for generating pattern of a given length
def seq_pattern(file1,file2,num):
    df1 = file1
    df1.columns = ['Seq']
    df2 = file2
    df2.columns = ['Name']
    cc = []
    dd = []
    ee = []
    for i in range(len(df1)):
        for j in range(len(df1['Seq'][i])):
            xx = df1['Seq'][i][j:j+num]
            if len(xx) == num:
                cc.append(df2['Name'][i])
                dd.append('Pattern_'+str(j+1))
                ee.append(xx)
    df3 = pd.concat([pd.DataFrame(cc),pd.DataFrame(dd),pd.DataFrame(ee)],axis=1)
    df3.columns= ['Seq ID','Pattern ID','Seq']
    return df3


# defining the function to process the blast output
def BLAST_processor(blast_result,name1,ml_results,thresh):
    name1.columns = [0]
    print(thresh)
    if os.stat(blast_result).st_size != 0:
        df1 = pd.read_csv(blast_result, sep="\t", names=['name','hit','identity','r1','r2','r3','r4','r5','r6','r7','r8','r9'])
        df__2 = name1
        df2 = pd.DataFrame()
        df2 = pd.concat([df2,df__2])
        df3 = ml_results
        cc = []
        for i in df2[0]:
            kk = i.replace('>','')
            if len(df1.loc[df1.name==kk])>0:
                df4 = df1[['name','hit']].loc[df1['name']==kk].reset_index(drop=True)
                if df4['hit'][0].split('_')[2]=='1':
                    cc.append(0.5)
                if df4['hit'][0].split('_')[2]=='0':
                    cc.append(-0.5)
            else:
                cc.append(0)
        df6 = pd.DataFrame()
        df6['Seq ID'] = [i.replace('>','') for i in df2.iloc[:,0]]
        df6['ML Score'] = df3['ML Score']
        df6['BLAST Score'] = cc
        df6['Hybrid Score'] = df6['ML Score']+df6['BLAST Score']
        df6['Prediction'] = ['ABPs' if df6['Hybrid Score'][i]>thresh else 'non-ABPs' for i in range(0,len(df6))]
    else:
        df__2 = name1
        df3 = ml_results
        df2 = pd.DataFrame()
        df2 = df2 = pd.concat([df2,df__2])
        ss = []
        vv = []
        for j in df2[0]:
            ss.append(j.replace('>',''))
            vv.append(0)
        df6 = pd.DataFrame()
        df6['Seq ID'] = ss
        df6['ML Score'] = df3['ML Score']
        df6['BLAST Score'] = vv
        df6['Hybrid Score'] = df6['ML Score']+df6['BLAST Score']
        df6['Prediction'] = ['ABPs' if df6['Hybrid Score'][i]>thresh else 'non-ABPs' for i in range(0,len(df6))]
    return df6

def BLAST_search(blast_result,name1):
    name1.columns = [0]
    if os.stat(blast_result).st_size != 0:
        df1 = pd.read_csv(blast_result, sep="\t", names=['name','hit','identity','r1','r2','r3','r4','r5','r6','r7','r8','r9'])
        df__2 = name1
        df2 = pd.DataFrame()
        df2 = pd.concat([df2,df__2])
        cc = []
        for i in df2[0]:
            kk = i.replace('>','')
            if len(df1.loc[df1.name==kk])>0:
                df4 = df1[['name','hit']].loc[df1['name']==kk].reset_index(drop=True)
                if df4['hit'][0].split('_')[2]=='1':
                    cc.append(1)
                if df4['hit'][0].split('_')[2]=='0':
                    cc.append(0)
            else:
                cc.append(0)
        df6 = pd.DataFrame()
        df6['Seq ID'] = [i.replace('>','') for i in df2.iloc[:,0]]
        df6['BLAST Score'] = cc
        df6['Prediction'] = ['ABPs' if df6['BLAST Score'][i]>0.5 else 'non-ABPs' for i in range(0,len(df6))]
    else:
        df__2 = name1
        df2 = pd.DataFrame()
        df2 = df2 = pd.concat([df2,df__2])
        ss = []
        vv = []
        for j in df2[0]:
            ss.append(j.replace('>',''))
            vv.append(0)
        df6 = pd.DataFrame()
        df6['Seq ID'] = ss
        df6['BLAST Score'] = vv
        df6['Prediction'] = ['ABPs' if df6['BLAST Score'][i]>0.5 else 'non-ABPs' for i in range(0,len(df6))]
    return df6


def MERCI_Processor(merci_file,merci_processed,name):
    hh =[]
    jj = []
    kk = []
    qq = []
    filename = merci_file
    df = pd.DataFrame(name)
    zz = list(df[0])
    check = '>'
    with open(filename) as f:
        l = []
        for line in f:
            if not len(line.strip()) == 0 :
                l.append(line)
            if 'COVERAGE' in line:
                for item in l:
                    if item.lower().startswith(check.lower()):
                        hh.append(item)
                l = []
    if hh == []:
        ff = [w.replace('>', '') for w in zz]
        for a in ff:
            jj.append(a)
            qq.append(np.array(['0']))
            kk.append('non-ABPs')
    else:
        ff = [w.replace('\n', '') for w in hh]
        ee = [w.replace('>', '') for w in ff]
        rr = [w.replace('>', '') for w in zz]
        ff = ee + rr
        oo = np.unique(ff)
        df1 = pd.DataFrame(list(map(lambda x:x.strip(),l))[1:])
        df1.columns = ['Name']
        df1['Name'] = df1['Name'].str.strip('(')
        df1[['Seq','Hits']] = df1.Name.str.split("(",expand=True)
        df2 = df1[['Seq','Hits']]
        df2.replace(to_replace=r"\)", value='', regex=True, inplace=True)
        df2.replace(to_replace=r'motifs match', value='', regex=True, inplace=True)
        df2.replace(to_replace=r' $', value='', regex=True,inplace=True)
        total_hit = int(df2.loc[len(df2)-1]['Seq'].split()[0])
        for j in oo:
            if j in df2.Seq.values:
                jj.append(j)
                qq.append(df2.loc[df2.Seq == j]['Hits'].values)
                kk.append('ABPs')
            else:
                jj.append(j)
                qq.append(np.array(['0']))
                kk.append('non-ABPs')
    df3 = pd.concat([pd.DataFrame(jj),pd.DataFrame(qq),pd.DataFrame(kk)], axis=1)
    df3.columns = ['Name','Hits','Prediction']
    df3.to_csv(merci_processed,index=None)

def Merci_after_processing(merci_processed,final_merci):
    df5 = pd.read_csv(merci_processed)
    df5 = df5[['Name','Hits']]
    df5.columns = ['Subject','Hits']
    kk = []
    for i in range(0,len(df5)):
        if df5['Hits'][i] > 0:
            kk.append(0.5)
        else:
            kk.append(0)
    df5["MERCI Score"] = kk
    df5 = df5[['Subject','MERCI Score']]
    df5.to_csv(final_merci, index=None)

def hybrid(ML_output,name1,merci_output,threshold,final_output):
    df6_2 = pd.read_csv(ML_output,header=None)
    df6_1 = pd.DataFrame(name1)
    df5 = pd.read_csv(merci_output)
    df4 = pd.read_csv(blast_output)
    df6 = pd.concat([df6_1,df6_2],axis=1)
    df6.columns = ['Subject','ML Score']
    df6['Subject'] = df6['Subject'].str.replace('>','')
    df7 = pd.merge(df6,df5, how='outer',on='Subject')
    df8 = pd.merge(df7,df4, how='outer',on='Subject')
    df8.fillna(0, inplace=True)
    df8['Hybrid Score'] = df8.sum(axis=1)
    df8 = df8.round(3)
    ee = []
    for i in range(0,len(df8)):
        if df8['Hybrid Score'][i] > float(threshold):
            ee.append('Toxin')
        else:
            ee.append('Non-Toxin')
    df8['Prediction'] = ee
    df8.to_csv(final_output, index=None)

# defining the function to read and implement the models
def model_run(file1,file2):
    a = []
    data_test = file1
    clf = pickle.load(open(file2,'rb'))
    y_p_score1=clf.predict_proba(data_test)
    y_p_s1=y_p_score1.tolist()
    a.extend(y_p_s1)
    df = pd.DataFrame(a)
    df1 = df.iloc[:,-1].round(2)
    df2 = pd.DataFrame(df1)
    df2.columns = ['ML Score']
    return df2
    
if __name__ == "__main__":
    print('############################################################################################')
    print('# This program AntiBP3 is developed for predicting and scanning          #')
    print('# Antibacterial peptides, developed by Prof G. P. S. Raghava group.               #')
    print('# Please cite: AntiBP3; available at https://webs.iiitd.edu.in/raghava/antibp3/  #')
    print('############################################################################################')

    parser = argparse.ArgumentParser(description='Please provide following arguments')

## Read Arguments from command
    parser.add_argument("-i", "--input", type=str, required=True, help="Input: protein or peptide sequence(s) in FASTA format or single sequence per line in single letter code")
    parser.add_argument("-o", "--output",type=str, help="Output: File for saving results by default outfile.csv")
    parser.add_argument("-s", "--source",type=int, choices = [1,2,3], help="Source: 1:GP ABPs, 2:GN ABPs, 3:GV ABPs by default 1")
    parser.add_argument("-j", "--job",type=int, choices = [1,2,3,4,5], help="Job Type: 1:Predict, 2:Design, 3:BLAST Search 4:Motif Scan, 5:Protein Scan ; by default 1")
    parser.add_argument("-t","--threshold", type=float, help="Threshold: Value between 0 to 1 by default 0.5 for GP ABPs, 0.45 for GN ABPs and 0.51 for GV ABPs")
    parser.add_argument("-e","--eval", type=float, help="E-value for Blast search (Blast Search only), by default 0.01 for GP ABPs, 0.01 for GN ABPs and 0.001 for GV ABPs")
    parser.add_argument("-w","--winleng", type=int, choices =range(8, 21), help="Window Length: 8 to 20 (scan mode only), by default 8")
    parser.add_argument("-d","--display", type=int, choices = [1,2], help="Display: 1:ABPs only, 2: All peptides, by default 1")
    parser.add_argument("-wd", "--working",type=str, help="Working Directory: Location for writing results")
    args = parser.parse_args()
    
    # Parameter initialization or assigning variable for command level arguments
    
    Sequence= args.input        # Input variable 
     
    # Output file 
     
    if args.output == None:
        result_filename= "outfile.csv" 
    else:
        result_filename = args.output
    # Source
    if args.source == None:
        source = int(1)
    else:
        source = args.source
    #Category
    if source == 1:
        category = 'GP ABPs'
    elif source == 2:
        category = 'GN ABPs'
    else:
        category = 'GV ABPs'
    
    # Threshold 
    if source == 1:
        if args.threshold is None:
            threshold = 0.5
        else:
            threshold = float(args.threshold)
    elif source == 2:
        if args.threshold is None:
            threshold = 0.45
        else:
            threshold = float(args.threshold)
    else:
        if args.threshold is None:
            threshold = 0.51
        else:
            threshold = float(args.threshold)
    # Job Type 
    if args.job == None:
        Job = int(1)
    else:
        Job = int(args.job)
    
    if source == 1:
        if args.eval is None:
            eval = 0.01
        else:
            eval = float(args.eval)
    elif source == 2:
        if args.eval is None:
            eval = 0.01
        else:
            eval = float(args.eval)
    else:
        if args.eval is None:
            eval = 0.001
        else:
            eval = float(args.eval)
    
    # Display
    if args.display == None:
        dplay = int(2)
    else:
        dplay = int(args.display)

    # Window Length 
    if args.winleng == None:
        Win_len = int(8)
    else:
        Win_len = int(args.winleng)

    # Working directory
    if args.working == None:
            wd = '.'
    else:
            wd = args.working


#####################################BLAST Path############################################
    if os.path.exists('envfile'):
        with open('envfile', 'r') as file:
            data = file.readlines()
        output = []
        for line in data:
            if not "#" in line:
                output.append(line)
        print(output)
        if len(output)==2:
            paths = []
            for i in range (0,len(output)):
                paths.append(output[i].split(':')[1].replace('\n',''))
            print(paths)
            blastp = paths[0]
            blastdb = paths[1]
            
        else:
            print("####################################################################################")
            print("Error: Please provide paths for BLAST, and required files", file=sys.stderr)
            print("####################################################################################")
            sys.exit()

    else:
        print("####################################################################################")
        print("Error: Please provide the '{}', which comprises paths for BLAST".format('envfile'), file=sys.stderr)
        print("####################################################################################")
        sys.exit()
#######################################################################################

    
    ##################################   BLAST Path    ############################################
    nf_path = os.path.dirname(__file__)
    blastdb1 = blastdb + "/grampos_db"
    blastdb2 = blastdb + "/gramneg_db"
    blastdb3 = blastdb + "/gramvariable_db"
    merci = nf_path + "/MERCI/MERCI_motif_locator.pl"
    merci_motif1 = nf_path + "/motif/GP_motif"
    merci_motif2 = nf_path + "/motif/GN_motif"
    merci_motif3 = nf_path + "/motif/GV_motif"
    ###########################################################################################
    
    if Job==3:
        print("\n");
        print('##############################################################################')
        print('Summary of Parameters:')
        print('Input File: ',Sequence,'; Threshold: ', threshold,'; Job Type: ', Job,';Category: ',category)
        print('Output File: ',result_filename,'; E-value ',eval,'; Display: ',dplay)
        print('##############################################################################')
    else:
        print("\n");
        print('##############################################################################')
        print('Summary of Parameters:')
        print('Input File: ',Sequence,'; Threshold: ', threshold,'; Job Type: ',Job)
        print('Output File: ',result_filename,'; Display: ',dplay,';Category: ',category)
        print('# ############################################################################')
        
    #======================= First module : Prediction  =====================#
    if Job == 1:
        print('\n======= Thanks for using Predict module of AntiBP3. Your results will be stored in file :',result_filename,' =====\n')
        df_2,dfseq = readseq(Sequence)
        print(df_2)
        df1 = lenchk(dfseq)
        X = feature_gen(df1, 8, wd)
        if source == 1:
            mlres = model_run(X, nf_path + '/model/modelRF_GP_aabNC.pkl')
        elif source == 2:
            mlres = model_run(X, nf_path + '/model/modelET_GN_aabNC.pkl')
        else:
            mlres = model_run(X, nf_path + '/model/modelSVC_GV_aabNC.pkl')
        filename = wd + "/" + str(uuid.uuid4())
        df11 = pd.concat([df_2,df1],axis=1)
        df11.to_csv(filename,index=None,header=False,sep="\n")
        mlres = mlres.round(3)
        df44 = mlres
        df44['Seq ID'] = [i.replace('>','') for i in df_2.iloc[:,0]]
        df44['Sequence'] = df1.Seq
        df44['Prediction'] = ['ABPs' if df44['ML Score'][i]>threshold else 'non-ABPs' for i in range(0,len(df44))]
        df44 = df44[['Seq ID','Sequence','ML Score','Prediction']]
        df44.to_csv(result_filename, index=None)
        print("\n=========Process Completed. Have a great day ahead.=============\n") 
    
        #===================== Second module : Design ======================#
    elif Job == 2:
        print('\n======= Thanks for using Design module of AntiBP3. Your results will be stored in file :',result_filename,' =====\n')
        print('==== Designing Peptides: Processing sequences please wait ...')
        df_2,dfseq = readseq(Sequence)
        df1 = lenchk(dfseq)
        df_1 = mutants(df1,df_2)
        dfseq = df_1[['Seq']]
        X = feature_gen(dfseq, 8, wd)
        if source == 1:
            mlres = model_run(X, nf_path + '/model/modelRF_GP_aabNC.pkl')
        elif source == 2:
            mlres = model_run(X, nf_path + '/model/modelET_GN_aabNC.pkl')
        else:
            mlres = model_run(X, nf_path + '/model/modelSVC_GV_aabNC.pkl')
        filename = wd + "/" + str(uuid.uuid4())
        df_1['Mutant'] = ['>'+df_1['Mutant ID'][i] for i in range(len(df_1))]
        df11 = df_1[['Mutant','Seq']] 
        df11.to_csv(filename,index=None,header=False,sep="\n")
        mlres = mlres.round(3)
        df44 = mlres
        df44['Mutant ID'] = [i.replace('>','') for i in df_1['Mutant']]
        #['_'.join(df44['Seq ID'][i].split('_')[:-1]) for i in range(len(df44))]
        df44['Seq ID'] = [i.replace('>','') for i in df_1['Seq ID']]
        df44['Sequence'] = df_1.Seq
        df44['Prediction'] = ['ABPs' if df44['ML Score'][i]>threshold else 'non-ABPs' for i in range(0,len(df44))]
        df44 = df44[['Seq ID','Mutant ID','Sequence','ML Score','Prediction']]
        if dplay == 1:
            df44 = df44.loc[df44.Prediction=="ABPs"]
        else:
            df44 = df44
        df44 = round(df44,3)
        df44.to_csv(result_filename, index=None)
        print("\n=========Process Completed. Have a great day ahead.=============\n")
        
    #=============== Third module : Blast Search ==================#
    if Job == 3:
        print('\n======= Thanks for using Blast scan module of AntiBP3. Your results will be stored in file :',result_filename,' =====\n')
        df_2,dfseq = readseq(Sequence)
        df1 = lenchk(dfseq)
        filename = wd + '/' +str(uuid.uuid4())
        df11 = pd.concat([df_2,df1],axis=1)
        df11.to_csv(filename,index=None,header=False,sep="\n")
        if source == 1:
            os.system(blastp + " -task blastp -db " + blastdb1 + " -query " + filename + " -out " + wd + "/RES_1_6_6.out -outfmt 6 -evalue " + str(eval))
        elif source == 2:
            os.system(blastp + " -task blastp -db " + blastdb2 + " -query " + filename + " -out " + wd + "/RES_1_6_6.out -outfmt 6 -evalue " + str(eval))
        else:
            os.system(blastp + " -task blastp -db " + blastdb3 + " -query " + filename + " -out " + wd + "/RES_1_6_6.out -outfmt 6 -evalue " + str(eval))
        df44 = BLAST_search(wd + '/RES_1_6_6.out',df_2)
        df44['Sequence'] = df1.Seq
        df44 = df44[['Seq ID','Sequence','BLAST Score','Prediction']]
        if dplay == 1:
            df44 = df44.loc[df44.Prediction=="ABPs"]
        else:
            df44 = df44
        df44 = round(df44,3)
        df44.to_csv(result_filename, index=None)
        print("\n=========Process Completed. Have a great day ahead.=============\n") 
    
    #=============== Fourth module : Motif Scan ==================#
    if Job == 4:
        print('\n======= Thanks for using Motif Scan module of AntiBP3. Your results will be stored in file :',result_filename,' =====\n')
        df_2,dfseq = readseq(Sequence)
        df1 = lenchk(dfseq)
        filename = wd + '/' + str(uuid.uuid4())
        df11 = pd.concat([df_2,df1],axis=1)
        df11.to_csv(filename,index=None,header=False,sep="\n")
        if source == 1:
            os.system("perl " + merci + " -p " + filename + " -i " + merci_motif1 + " -o " + wd + "/merci.out")
        elif source == 2:
            os.system("perl " + merci + " -p " + filename + " -i " + merci_motif2 + " -o " + wd + "/merci.out")
        else:
            os.system("perl " + merci + " -p " + filename + " -i " + merci_motif3 + " -o " + wd + "/merci.out")    
        MERCI_Processor(wd + "/merci.out", wd + "/merci_processed_out", df_2)
        df44 = pd.read_csv(wd + "/merci_processed_out")
        if dplay == 1:
            df44 = df44.loc[df44.Prediction=="ABPs"]
        else:
            df44 = df44
        df44 = round(df44,3)
        df44.to_csv(result_filename, index=None)
        print("\n=========Process Completed. Have a great day ahead.=============\n")

    #=============== Fifth module : Protein Scan ==================#
    if Job == 5:
        print('\n======= Thanks for using Protein Scan module of AntiBP3. Your results will be stored in file :',result_filename,' =====\n')
        print('==== Scanning Peptides: Processing sequences please wait ...')
        df_2,dfseq = readseq(Sequence)
        df_1 = seq_pattern(dfseq,df_2,Win_len)
        dfseq = df_1[['Seq']]
        X = feature_gen(dfseq, 8, wd)
        if source == 1:
            mlres = model_run(X, nf_path + '/model/modelRF_GP_aabNC.pkl')
        elif source == 2:
            mlres = model_run(X, nf_path + '/model/modelET_GN_aabNC.pkl')
        else:
            mlres = model_run(X, nf_path + '/model/modelSVC_GV_aabNC.pkl')
        filename = str(uuid.uuid4())
        df_1['Pattern'] = ['>'+df_1['Pattern ID'][i] for i in range(len(df_1))]
        df11 = df_1[['Pattern','Seq']]
        df11.to_csv(filename,index=None,header=False,sep="\n")
        mlres = mlres.round(3)
        df44 = mlres
        df44['Pattern ID'] = [i.replace('>','') for i in df_1['Pattern']]
        df44['Seq ID'] = [i.replace('>','') for i in df_1['Seq ID']]
        df44['Sequence'] = df_1.Seq
        df44['Prediction'] = ['ABPs' if df44['ML Score'][i]>threshold else 'non-ABPs' for i in range(0,len(df44))]
        df44 = df44[['Seq ID','Pattern ID','Sequence','ML Score','Prediction']]
        if dplay == 1:
            df44 = df44.loc[df44.Prediction=="ABPs"]
        else:
            df44 = df44
        df44 = round(df44,3)
        df44.to_csv(result_filename, index=None)
        print("\n=========Process Completed. Have a great day ahead ahead.=============\n")

