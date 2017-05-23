#!/usr/bin/python
#
####################################################
#procedure para_input :
# Input the parameter for further analysis

def para_input(hd_cutoff,read_comb,truncate,singleton,infile1,infile2,g_start,g_end) :

 infile_read = open('parameter.txt','r')

 for ln in infile_read.readlines() :
  array = ln[0:len(ln)-1].split()
  if array[0] == "hd_cutoff" : # The similarity cutoff for hamming distance between two sequences
   hd_cutoff = int(array[1])
  if array[0] == "read_comb" : # The label to indicate whether the read1 and read2 should be combined
                               # 0. combine read1 and read2, 1. read1, 2. read2
   read_comb = int(array[1])
   if read_comb == 0 or read_comb == 1 :
    read_reverse = 1  # The labe to indicate whether the reads need to be reversely complemented
   else :
    read_reverse = 0
  if array[0] == "truncate" : # The truncation length of reads. If truncate = 0, just use the original length of reads
   truncate = int(array[1])
  if array[0] == "singleton" :
   singleton = int(array[1]) # The labe to indicate whether the analysis will be made for singleton reads
                             # 0 indicates do not make the singleton analysis
                             # 1 indicates make the singleton analysis
  if array[0] == "read1" :
   infile1 = array[1]
  if array[0] == "read2" :
   infile2 = array[1]
  if array[0] == "g_start" :
   g_start = int(array[1])
  if array[0] == "g_end" :
   g_end = int(array[1])
 return (hd_cutoff,read_comb,truncate,singleton,infile1,infile2,g_start,g_end)
##############################################################
#procedure hamming :
# calculate the hamming distance

def hamming(string1,string2) :
 number = 0
 if len(string1) != len(string2) :
  number = 1000
  return number
 else :
  for i in range(0,len(string1)) :
   if string1[i:i+1] != string2[i:i+1] :
    number = number + 1                                                                                           
  return number
#####################################################
#procedure levenshtein :
def levenshtein(a,b):
    "Calculates the Levenshtein distance between a and b."
    n, m = len(a), len(b)
    if n > m:
        # Make sure n <= m, to use O(min(n,m)) space
        a,b = b,a
        n,m = m,n

    current = range(n+1)
    for i in range(1,m+1):
        previous, current = current, [i]+[0]*n
        for j in range(1,n+1):
            add, delete = previous[j]+1, current[j-1]+1
            change = previous[j-1]
            if a[j-1] != b[i-1]:
                change = change + 1
            current[j] = min(add, delete, change)

    return current[n]
#####################################################
#procedure pairwise :

def pairwise(number,barcode,temp,pair_dist) :

# out_pair = open(barcode+".res",'w')
 iden = {}

 for i in range(1,number) :
  for j in range(i+1,number+1) :
   tag_1 = barcode + "_" + str(i) + "_" + str(j)
   tag_2 = barcode + "_" + str(j) + "_" + str(i)
   if i in iden.keys() :
    tag_3 = barcode + "_" + str(iden[i]) + "_" + str(j)
    pair_dist[tag_1] = pair_dist[tag_3]
    pair_dist[tag_2] = pair_dist[tag_3]
   else :
    pair_dist[tag_1] = levenshtein(temp[barcode+"_"+str(i)],temp[barcode+"_"+str(j)])
    pair_dist[tag_2] = pair_dist[tag_1]
#   out_pair.write(str(i)+"\t"+str(j)+"\t"+str(pair_dist[tag_1])+"\n")
   if int(pair_dist[tag_1]) == 0 :
    if j not in iden.keys() :
     iden[j] = i
 return (pair_dist)
#####################################################
#procedure cluster :                                                                                              
#def p_cluster(number,bar,temp,temp_score,hd_cutoff,file_out,singleton) :                                                                                                                                  
def p_cluster(number,bar,temp,temp_score,hd_cutoff,file_out,lst_out,singleton,pair_dist) :                                          
                                                                                                                  
# number : the number of raw reads with the identical barcode                                                     
# barcode : the barcode for clustering                                                                            
# temp : the temporary data for raw reads                                                                         
# temp_score : the temporary data for the score of raw reads                                                      
# hd_cutoff : the similarity distance cutoff for two sequencing, using hamming distance                           
# file_cluster : the file to record the clustering and build consensus results                                    
# file_singleton : the file to record the singleton cluster results                                               
 group = {}
 cons = {}
 value = {}
 ss = {}
 consensus = ""
 single_cons = ""
 sensus = 0
 clustering = "True"
                                                                                    
# Cluster
# start = time.clock()
 length = len(temp[bar+"_1"])            
 if number == 1 :                                                                 
   group[bar+"_"+str(number)]=bar+"_"+str(number)             
 else :   
   # First step : find the seed candicate with the most similar reads with distance < threshold
   # for two reads with similarity < hamming distance cutoff                      
   # cluster together, use the barcode(bar+"_"+str(i)) to label these two reads      
   for i in range(1,number) :
    for j in range(i+1,number+1) :
     clustering = "True"
     tag = bar+"_"+str(i)+"_"+str(j)
     tag_lr = bar+"_"+str(i)   #tag_lr : label reference
     tag_lq = bar+"_"+str(j)   #tag_lq : label query
     if int(pair_dist[tag]) <= hd_cutoff :
      for k in range(1,i) :
       tag_lk = bar+"_"+str(k)
       tag_pi = bar+"_"+str(j)+"_"+str(k)
       if group[tag_lk] == group[tag_lr] :
        if pair_dist[tag_pi] > hd_cutoff :
         clustering = "False"
         break
     else :
      clustering = "False"
     if clustering == "True" :
      if bar+"_"+str(i) in group :            
       if group[bar+"_"+str(i)] == bar+"_"+str(i) :                 
        group[bar+"_"+str(j)]=bar+"_"+str(i)                      
       else :                                                 
        group[bar+"_"+str(j)]=group[bar+"_"+str(i)]                
      else :                               
       group[bar+"_"+str(i)]=bar+"_"+str(i)                          
       group[bar+"_"+str(j)]=bar+"_"+str(i)                       
     else :
      if bar+"_"+str(i) not in group :
       group[bar+"_"+str(i)]=bar+"_"+str(i)
      if bar+"_"+str(j) not in group :
       group[bar+"_"+str(j)]=bar+"_"+str(j)
# elapsed = (time.clock() - start )
# print "Cluster","\t","time","\t",elapsed                                                                        

# Build consensus                                                                                                 
# start = time.clock()
 for i in group.values() :                                                                                        
    if i not in value :                                                                                            
     value[i] = 1                                                                                                  
    else :                                                                                                         
     value[i] = value[i] + 1                                                                                       
 for i in value.keys() :
    length = len(temp[i])
    for k in range(0,length) :                                                                                   
     for j in range(1,number+1) :                                                                                  
      if group[bar+"_"+str(j)] == i :
       if temp[bar+"_"+str(j)][k] not in cons :                                                                
        cons[temp[bar+"_"+str(j)][k]]=1                                                                        
        ss[temp[bar+"_"+str(j)][k]] = ord(temp_score[bar+"_"+str(j)][k])-33                                
       else :                                                                                                      
        cons[temp[bar+"_"+str(j)][k:k+1]]=cons[temp[bar+"_"+str(j)][k:k+1]]+1                                      
        ss[temp[bar+"_"+str(j)][k:k+1]] = ss[temp[bar+"_"+str(j)][k:k+1]]+ord(temp_score[bar+"_"+str(j)][k:k+1])-33
     for l in cons.keys() :                                                                                        
      if cons[l] == max(cons.values()) and ss[l] == max(ss.values()) :                                             
#      if cons[l] == max(cons.values()) :
       consensus = consensus + l                                                                                   
       if singleton == 1 :                                                                                         
        sensus = sensus - cons[l]                                                                                  
        single_cons = single_cons + chr(sensus+33)                                                                 
       break                                                                                                       
     cons.clear()                                                                                                  
     ss.clear()                                                                                                    
     sensus = 0                                                                                                    

    for j in range(1,number+1) :
     if group[bar+"_"+str(j)] == i :
      lst_out.write(">"+i+"_"+str(value[i])+"\t"+str(j)+"\n")
                                                                                                                  
    if singleton == 0 :                                                                                            
     file_out.write(">"+i+"_"+str(value[i]))                                                                       
     file_out.write("\n")                                                                                          
     file_out.write(consensus)                                                                                     
     file_out.write("\n")                                                                                          
    if singleton == 1 :                                                                                            
     file_out.write(">"+i+"_"+str(value[i]))                                                                       
     file_out.write("\n")                                                                                          
     file_out.write(consensus)                                                                                     
     file_out.write("\n")                                                                                          
     file_out.write(single_cons)                                                                                   
     file_out.write("\n")                                                                                          
                                                                                                                  
    consensus = ""                                                                                                 
    single_cons = ""                                                                                               
    sensus = 0                                                                                                     
 group.clear()                                                                                                    
 cons.clear()                                                                                                     
 value.clear()                                                                                                    
 ss.clear()
# elapsed = (time.clock() - start )
# print "Consensus","\t","time","\t",elapsed
######################################################                                                            

import time

hd_cutoff = 10
read_comb = 0
truncate = 0
singleton = 0
infile1 = ""
infile2 = ""
g_start = 0
g_end = 0

temp = {}
temp_score = {}
pair_dist = {}

(hd_cutoff,read_comb,truncate,singleton,infile1,infile2,g_start,g_end) = para_input(hd_cutoff,read_comb,truncate,singleton,infile1,infile2,g_start,g_end)

single = open('singleton.fasta','w')

for i in range(g_start,g_end+1) :
# print i
 number = 0
 count = 0
 count_2 = 0
 infile = open('group_'+str(i)+'.txt','r')
 file_cluster = open('group_'+str(i)+'_cluster.txt','w')
 cluster_lst = open('group_'+str(i)+'_cluster.lst','w')
 for ln in infile.readlines() :
  if ln[0] == ">" :
   if count == 1 :
    single.write(">" + barcode + "_10000_1" + "\n")
    single.write(temp[barcode+"_1"]+"\n")
   if count > 1 :
    if singleton == 0 :
     pair_dist = pairwise(number,barcode,temp,pair_dist)
     p_cluster(number,barcode,temp,temp_score,hd_cutoff,file_cluster,cluster_lst,singleton,pair_dist) #Call the process to cluster the barcode
     temp.clear()
     temp_score.clear()
     pair_dist.clear()
   count = 0
   number = 0
   barcode = ln.strip()[1::]
   temp.clear()
  else :
   array = ln.strip().split("|")
   count = count + 1
   number = number + 1
   temp[barcode+"_"+str(number)]=array[0]
   temp_score[barcode+"_"+str(number)] = array[1]

 if count == 1:
    single.write(">" + barcode + "_10000_1" + "\n")
    single.write(temp[barcode+"_1"]+"\n")
 if count >1 :
  pair_dist = pairwise(number,barcode,temp,pair_dist)
  p_cluster(number,barcode,temp,temp_score,hd_cutoff,file_cluster,cluster_lst,singleton,pair_dist)
  temp.clear()
  temp_score.clear()
  pair_dist.clear()
  count = 0
  number = 0
  
 infile.close()
 file_cluster.close()
