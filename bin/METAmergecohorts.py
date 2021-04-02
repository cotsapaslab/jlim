
import numpy as np
import scipy.stats as sta
import time
import os, sys
import gzip
import os
import shutil


mergenum=int(sys.argv[1])   # how many cohorts are you merging
assocfile=sys.argv[2] 
outfile=sys.argv[3]   # name of output file

Betas=[]
Vars=[]

for k in range(mergenum):
    Betas.append( np.loadtxt(sys.argv[2*k+4]) )
    Vars.append( np.loadtxt(sys.argv[2*k+5]) )

nperms=len(Betas[0])  # dimension of matrix, really nperm+1
nloci=len(Betas[0][0])

totW=np.zeros([nperms,nloci]) # total sum weight
metaB=np.zeros([nperms,nloci]) # total sum weight
weight=np.empty([nperms,nloci])

for k in range(mergenum):
    weight=1/Vars[k]
    totW+=weight
    metaB+=Betas[k]*weight
tstats=metaB/np.sqrt(totW)   #This is /totw for weight corection times sqrt(totw) for division of meta-standard deviation 

BETA=metaB[0]/totW[0]
#pvs=2-2*sta.norm.cdf(abs(tstats[0]))   # cdf function is not very good at low p values- underflow
pvs=2*sta.norm.sf( abs(tstats[0]) )   # sf- survival function, 1-cdf is way more precise, can go down to -logp of 280

outfilename = outfile +'.meta.mperm.dump.all.gz'
output = gzip.open(outfilename, 'wb')
try:
    for j in range(len(tstats)):
        output.write(str(j)+" "+" ".join([str(i)[:8] for i in tstats[j]])+"\n")
finally:
    output.close()

f=gzip.open(assocfile)
lines=f.readlines()
f.close()

CHR=lines[1].split(" ")[0]
nmiss=lines[1].split(" ")[5]
snp=[]
bp=[]
bimline=[]
for k in range(1,len(lines)):
    vec=lines[k].split(" ")
    snp.append(vec[1])
    bp.append(vec[2])
    bimline.append(vec[3])
    
outfilenameA = outfile +'.meta.assoc.linear.gz'
output = gzip.open(outfilenameA, 'wb')
# Separate entries with a single space
output.write("CHR SNP BP BIMLINE TEST NMISS BETA STAT P\n")  # BIMline is the corresponding line in bimbam file from chromosome
for loc in range(nloci):
    output.write(CHR+" "+snp[loc]+" "+bp[loc]+" "+bimline[loc]+" ADD "+nmiss+" "+str(BETA[loc])+" "+str(tstats[0][loc])+" "+str(pvs[loc])+"\n")
output.close()
