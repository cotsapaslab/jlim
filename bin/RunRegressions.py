# on clusters sometimes you need: export OPENBLAS_NUM_THREADS=1
 
import numpy as np
import time
from scipy.stats import t
import os, sys
#import statsmodels.api as sm
import gzip


genfile=sys.argv[1]
snpIDfile=sys.argv[2]  # ends in snps.gz
sampsfile=sys.argv[3]
phenfile=sys.argv[4]
covfile=sys.argv[5]
outfile=sys.argv[6]   # MrOS.1.135555.155555.rdi3p
covarlist=[int(k) for k in sys.argv[7].split(",")]
nperm=int(sys.argv[8])
CHR=sys.argv[9]
cohn=sys.argv[10]
 

positions=np.load(snpIDfile[:-7]+"merge"+cohn+".positions.npy") # This array holds the indexes of the SNPs that merged succesfully accross all cohorts
nloci=len(positions) 
flipped=np.arange(nloci)[positions<0]   # alleles coded in reverse
positions[flipped]=-positions[flipped]

f=open(phenfile)
lines=f.readlines()
f.close()
nsamples=len(lines)-1
sampleNames=[]
pheno=np.empty(nsamples,dtype=float)
for k in range(1,nsamples+1):
    pheno[k-1]=float(lines[k].split("\t")[-1][:-1])
    sampleNames.append(lines[k].split("\t")[0])
    
    
sampleIndeces=np.arange(nsamples)[np.invert(np.isnan(pheno))] #indeces of samples with phenotype in 
Phenotype=pheno[sampleIndeces]
nsamples=len(sampleIndeces)
print str(nsamples)+" samples" 

phenosamps=np.array(sampleNames)[sampleIndeces]

f=open(sampsfile)
lines=f.readlines()
f.close()
genosamps=[k.split(" ")[0] for k in lines]
genosamps=[k[:-1] if k[-1]=='\n' else k for k in genosamps]  # removes newline at end of sample name if it exists
try:
    genoSampIndeces=np.array([genosamps.index(k) for k in phenosamps])
except:
    print 'ERROR, sample in phenotype file missing from sample file. Check that format is identical. Sample file should have a single column'
    genoSampIndeces=np.array([genosamps.index(k) for k in phenosamps])

Loci=np.loadtxt(genfile)

if(covarlist[0]==0):
    CovMatrix=np.ones([2,nsamples])
    print "No covariates given"
else:
    Numcovs=len(covarlist)
    f=open(covfile)
    ncl=f.readlines()
    f.close()
    hh=ncl[0].split("\t")
    print 'Covariates'
    print [hh[k] for k in covarlist]
    CovMatrix=np.ones([Numcovs+2,nsamples])
    covsamps=[k.split("\t")[0] for k in ncl]
    try:
        covSampIndeces=np.array([covsamps.index(k) for k in phenosamps])
    except:
        print 'ERROR, sample in phenotype file missing from covariate file. Check that format is identical. Both should be tab separated'
        covSampIndeces=np.array([covsamps.index(k) for k in phenosamps])
    for i in range(nsamples):
        entry=ncl[covSampIndeces[i]].split("\t")
        for k in range(Numcovs):
            CovMatrix[k+2,i]=float(entry[covarlist[k]])

#print Counter([sampleNames[sampleIndeces[i]]==ids[i] for i in range(nsamples)])

# This version runs the regression only once per locus, then stores the (XTX)^(-1)XT matrix and the standard errors
# of the first coeficient. Then for each permutation it is faster to obtain parameters and t statistics.
# It also creates all the permutations once and stores them in newvecS.

st=time.time()

newvecS=np.empty([nperm,nsamples],dtype=int)
allV=np.arange(nsamples)
for k in range(nperm):   # this loop creates the random shuffeld ordering of the phenotypes
    np.random.shuffle(allV)
    newvecS[k]=np.array(allV)

BETA=[]
tSTAT=[]
pvs=[]    
BetaStat=np.empty([nperm+1,nloci])
varsStat=np.empty([nperm+1,nloci])
X=CovMatrix.T
kk=len(CovMatrix)
for loc in range(nloci):
    X[:,0]=Loci[positions[loc]][genoSampIndeces]
    xtxinv=np.linalg.inv(X.T.dot(X))
    xtxinvxt= xtxinv.dot(X.T)
    B=xtxinvxt.dot(Phenotype)
    BETA.append(B[0])
    stdcoef=np.sqrt(xtxinv[0,0])
    eo=( ((X.dot(B)-Phenotype).var())*nsamples/(nsamples-2.0) )**(0.5)  # estimate of error (SSE/n-2)
    tSTAT.append(B[0]/(stdcoef*eo))
    pvs.append(2*t.sf(abs(tSTAT[-1]),nsamples-kk))
    BetaStat[0,loc]=B[0]
    varsStat[0,loc]=(stdcoef*eo)**2
    for k in range(nperm):
        B=xtxinvxt.dot(Phenotype[newvecS[k]])
        e=( ((X.dot(B)-Phenotype[newvecS[k]]).var())*nsamples/(nsamples-2.0) )**(0.5)
        BetaStat[k+1,loc]=B[0]
        varsStat[k+1,loc]=(stdcoef*e)**2

BETA=np.array(BETA)
if(len(flipped)):
    BetaStat[:,flipped]=-BetaStat[:,flipped]
    BETA[flipped]=-BETA[flipped]

print "permutation statistics done. Time: "+str(time.time()-st)+"  Loci: "+str(nloci)

np.savetxt(outfile +'.betas.mperm.dump.all.gz',BetaStat,'%f8')
np.savetxt(outfile +'.vars.mperm.dump.all.gz',varsStat,'%f8')

gg=gzip.open(snpIDfile)
lin=gg.readlines()
gg.close()

snp=lin[0].split(",")
snp[-1]=snp[-1][:-1]
bp=lin[1].split(",")
bimline=lin[2].split(",")
bp[-1]=bp[-1][:-1]
bimline[-1]=bimline[-1][:-1]
nmiss=str(nsamples)

SE=np.sqrt(varsStat[0])

outfilenameA = outfile +'.assoc.linear.gz'
output = gzip.open(outfilenameA, 'wb')
#separate entries with a single space
output.write("CHR SNP BP BIMLINE TEST NMISS BETA SE STAT P\n")  # changed AIfor BIMline, the corresponding line in bimbam file from chromosome
for loc in range(nloci):
    output.write(CHR+" "+snp[positions[loc]]+" "+bp[positions[loc]]+" "+bimline[positions[loc]]+" ADD "+nmiss+" "+str(BETA[loc])+" "+str(SE[loc])+" "+str(tSTAT[loc])+" "+str(pvs[loc])+"\n")
output.close()
