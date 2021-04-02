
# Removes BOTH copies of repeated (multiallelic) SNPS which match at a position, even if rsid is different. 
# This version sorts rsids by the order in which they appear in original snp file in first cohort

import numpy as np
import os, sys
import gzip
from collections import Counter
from itertools import combinations

mergenum=int(sys.argv[1])   # how many cohorts are you merging
chro=sys.argv[2]  # chromosome
outfile=sys.argv[3]   # name of output file
bimbamfiles=[]
snpfiles=[]

for k in range(mergenum):
    bimbamfiles.append( sys.argv[2*k+4] )
    snpfiles.append( sys.argv[2*k+5] )

rsids=[]
posits=[]
alleles=[]

rsids2=[]
posits2=[]  # positions after removing duplicates (triallelic snps)

for j in range(mergenum):
    f=gzip.open(snpfiles[j])
    lines=f.readlines()
    f.close()
    rsids.append(lines[0].split(","))
    posits.append(np.array(lines[1].split(","),dtype=int))
    alleles.append(lines[3].split(","))
    rsids[-1][-1]=rsids[-1][-1][:-1]  # removes newline at end
    alleles[-1][-1]=alleles[-1][-1][:-1]
    assert len(posits[-1])==len(rsids[-1])
    cc=Counter(posits[-1])
    lc=list(cc)
    repeats=[]
    for i in range(len(lc)):
        if cc[lc[i]] !=1:
            repeats.append(lc[i])
    allpos=np.arange(len(posits[-1]))
    for k in repeats:
        allpos[np.arange(len(posits[-1]))[np.array(posits[-1])==k]]=-1
    if( len(repeats) ) :    
        allpos=np.unique(allpos)[1:]
    else:
        allpos=np.unique(allpos)
    rsids2.append(np.array(rsids[-1])[allpos])  
    rsidL=np.array(rsids2[-1])
    cr=Counter(rsidL)
    lc=list(cr)
    repeatsRsid=[]     # Checks if there are repeated rsid's that have differet bp coordinates
    for i in range(len(lc)):
        if cr[lc[i]] !=1:
            repeatsRsid.append(lc[i])
    for k in repeatsRsid:    # Changes the names of the repeated rsid's
        kpos=np.arange(len(rsids[-1]))[np.array(rsids[-1])==k]
        kpos2=np.arange(len(rsidL))[rsidL==k]
        for p in range(len(kpos)):
            rsids[-1][kpos[p]]=k+"v"+str(p)     # Makes the same name change in both rsids and rsids2
            rsids2[-1][kpos2[p]]=k+"v"+str(p)



comrsids=rsids2[0]
for j in range(mergenum-1):
    comrsids= list( set(comrsids) & set( rsids2[j+1] ))
    
totpos=len(comrsids) 
mergeindexR=np.empty([mergenum,totpos],dtype= int)   # the indeces of the positions from each snpfile (cohort) that will make it to the merged file- on posits[j]

mergeindexR[0]=np.array([rsids[0].index(rs) for rs in comrsids])
mergeindexR[0].sort()
sortcomrsids=[rsids[0][k] for k in mergeindexR[0] ]

for j in range(1,mergenum):
    mergeindexR[j]=np.array([rsids[j].index(rs) for rs in sortcomrsids])

mismatch=[]   #This will will contain indeces (in mergeindexR[j]) of mismatched alleles

refalleles=np.array(alleles[0])[ mergeindexR[0] ]
revs=0
for j in range(mergenum-1):
    for k in range(totpos):
        if not (refalleles[k] == alleles[j+1][ mergeindexR[j+1][k] ]):
            pair=alleles[j+1][ mergeindexR[j+1][k] ].split(".")
            mismatch.append(k)
            if( (pair[1]+"."+pair[0]) == refalleles[k]):    # IF true, this SNP is switched, ref<-->alternative
                print "reversed allele "+ refalleles[k]+ " , " + alleles[j+1][ mergeindexR[j+1][k]] + " cohort " +str(j+1) +" somefileinlocus "+ bimbamfiles[0]
                revs+=1


mismatch=np.unique(mismatch)
print "Variants matching in all cohorts: " + str(totpos)+ ", mismatched alleles: "+str(len(mismatch))+", reversed: " +str(revs)

totpos=totpos-len(mismatch)
mergeindexRr=np.empty([mergenum,totpos],dtype= int) 
for j in range(mergenum):
    mergeindexRr[j]=np.delete(mergeindexR[j],mismatch)
     

mLoci=[]
c=0
ee=[]
freqs=[]
for j in range(mergenum):
    mLoci.append( np.loadtxt(bimbamfiles[j] ) )
    le=len(mLoci[-1][0])
    ee.extend( ['samp'+str(i+1) for i in range(c,c+le) ] )
    c+=le
    freqs.append(np.array([j.mean()/2 for j in mLoci[-1]]))
    
output = gzip.open(outfile, 'wb')
output.write("clone\tChr\tStart\tEnd\t"+"\t".join(ee)+"\n")

for pos in range(totpos):
    e=[]
    physpo=np.array(posits[0],dtype=int)
    for j in range(mergenum):
        e.extend ([str(i) for i in mLoci[j][ mergeindexRr[j][pos] ,:]])
    output.write(rsids[0][mergeindexRr[0][pos]]+"\t"+chro+"\t"+str(physpo[mergeindexRr[0]][pos])+"\t"+str(physpo[mergeindexRr[0]][pos])+"\t"+"\t".join(e)+"\n")

output.close()

for j in range(mergenum):
    np.save(snpfiles[j][:-7]+"merge"+str(mergenum)+".positions",mergeindexRr[j])

kk=combinations(np.arange(mergenum),2)
for i in kk:
    dist=abs(freqs[i[0]][mergeindexRr[i[0]]]- freqs[i[1]][mergeindexRr[i[1]]]).mean()
    com=len( set( rsids2[i[0]]) & set(  rsids2[i[1]] ))
    print "Total rsids in cohort "+str(i[0]+1)+" : " + str(len( rsids2[i[0]]) ) +", and cohort "+str(i[1]+1)+ " : " + str(len( rsids2[i[0]]) ) + ". Common rsids: "+str(com)+ ", L1-distance " +str(dist)
print "first file " +bimbamfiles[0]
