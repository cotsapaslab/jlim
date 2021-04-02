
import numpy as np
import time
import os, sys
import gzip
import os
import shutil

bimbamfile=sys.argv[1]
mapfile=sys.argv[2]
peaksfile=sys.argv[3]
sampsfile=sys.argv[4]
outfile=sys.argv[5]

bimbamseparator=sys.argv[6]
mapseparator=sys.argv[7]
infopos=int(sys.argv[8]) # in map file (also called info file) index where chromosome position lies

try:
    mafT=float(sys.argv[9])  # Optional input MAF
except:
    mafT=0.001

try:
    interval=int(sys.argv[10])  # Optional input interval in bps
except:
    interval=200000

f=open(sampsfile)
lines=f.readlines()
f.close()
nsamps=len(lines)

f=gzip.open(mapfile)
lines=f.readlines()
f.close()

try:
    int(lines[0].split(mapseparator)[infopos])
except:
    assert 1==0,"ERROR mapfile and bimbam should not have a header, mapPosition should denote column in mapfile which holds genomic position, where the first column is column 0. Also check the mapfile separator"

vv=['0' for i in range(5)]
vv[infopos]='9000000000'
lines.append(mapseparator.join(vv)) # Cap introduced to stop loop below at the end of reading the file, in case that the last peak contains the last snp.

f=open(peaksfile)
pp=f.readlines()
f.close()
p=pp[0].split(",")
peaks=np.array(p,dtype=int)
peaks+=-interval/2
peaks.sort()

LP=len(peaks)



linepos=np.empty(len(peaks),dtype=int)    # will store the index of the line containing the first snp in the window that starts at position peak[k]
coun=0
pp=[]
print len(lines)
for i in range(len(lines)):
    pos=int(lines[i].split(mapseparator)[infopos])
    pp.append(pos)
    if pos>=peaks[coun]:
        linepos[coun]=i
        coun+=1
        if coun==LP:
            if len(pp)==1:  # introduced for a special case of single peak in the file causing a bug below
                pos=int(lines[i+1].split(mapseparator)[infopos])
                pp.append(pos)
            break


assert pp[0]!=pp[1],"ERROR  mapPosition should denote column in mapfile which holds genomic position, where the first column is column 0"
# This should happen if you pick the mapPosition correpsonding to the chromosome


snp=[]
abspos=[]   # Will have the positions of all the lines that carry SNPS within the windows in each peak
endpos=np.empty(len(linepos),dtype=int) # These will store the line number of the last snp that lies inside of each window
physpos=[]
for i in range(len(linepos)):
    c=0
    while(1):
        vec=lines[linepos[i]+c].split(mapseparator)   # 0 is snp ID, 1 is position, 2 is chromosome, (usually?)
        if ( (int(vec[infopos])>peaks[i]+interval)):  # if position is beyond interval we can go to next window
            endpos[i]=linepos[i]+c-1   # -1 to have last snp
            break
        snp.append(vec[0])
        abspos.append(linepos[i]+c)  
        physpos.append(vec[infopos])
        c+=1
abspos,indeces =np.unique(np.array(abspos),return_index=True)   # simplifies in case of overlapping peaks
physpos=np.array(physpos,dtype=int)[indeces]
# the unique function sorts and gets unique entries. with return index it also tells us the indeces of what it is returning (kind)
# of similar to np.argsort. This way physpos and abspos stay compatible

c=0
k=0
finsnp=[]      # Contains the SNP iDs
finpos=[]          #Cotains the line numbers corresponding to the final stored snps in the file
finphyspos=[]     # COntains physiscal position of snps
alleles=[]
lastpos=len(abspos)
arr=np.empty([1,nsamps],dtype=float)
Loci=np.empty([len(abspos),nsamps],dtype=float)
fp=gzip.open(bimbamfile)
for i, lineB in enumerate(fp):
    if(not i%100000):
        print i
    if i == abspos[c]:
        c+=1
        vec=lineB.split(bimbamseparator)
        arr[0]=np.array(vec[3:])    # only includes samples with phenotype. An error here can be caused by the bimbam file and sample file having different numbers of samples. These files should match exactly!
        mean=arr[0].mean()*0.5
        if mean>mafT and mean<(1-mafT):
            finsnp.append(vec[0])
            finpos.append(i)
            finphyspos.append(physpos[c-1])
            alleles.append(vec[1]+"."+vec[2])
            Loci[k,]=arr
            k+=1
        if(c==lastpos):
            break    #IF THIS IS true the last peak is inlcuded and we can now quit
finpos=np.array(finpos)
fp.close()

startindex=np.searchsorted(finpos,linepos-1)
endindex=np.searchsorted(finpos,endpos+1)

#These arrays tell us where each peak starts and ends in the arrays and matrices that exculsively have the final SNPS we are including

for p in range(len(peaks)):
    sliceMat=Loci[startindex[p]:endindex[p]]
    filename=outfile+"."+str(peaks[p])+"."+str(peaks[p]+interval)
    np.savetxt(filename+".data",sliceMat,fmt="%1.4f") # float, 4 decimal points
    f_in=open(filename+".data", 'rb')
    f_out=gzip.open(filename+".data.gz", 'wb')   # These lines gzip file
    shutil.copyfileobj(f_in, f_out)
    f_in.close()
    f_out.close()
    f=gzip.open(filename+".snps.gz",'wb')
    os.remove(filename+".data")
    f.write(",".join(finsnp[startindex[p]:endindex[p]])+"\n")
    f.write(",".join([str(s) for s in (finphyspos[startindex[p]:endindex[p]])])+"\n")
    f.write(",".join([str(s) for s in (finpos[startindex[p]:endindex[p]])])+"\n")
    f.write(",".join(alleles[startindex[p]:endindex[p]])+"\n")
    f.write("SNPID,Position in chromosome,corresponding line in complete BimBam file,reference-alternative")
    f.close()
