#reading in data funtion, basically want to read in different files line by line. Want to get rid of first line, not needed.
def readIn (path):
    input = open (path)
    lines =input.readlines()
    del lines[0]
    input.close()
    return lines
#Finding significant results accrossed datasets
#Lefse (at least) has a different data structure / critera than the rest of the
#the methods. Features with LDA scores >2.0 and adjusted P-values <.05 must be
# pulled. First read in the file, add ASV sequences to a list of significant resultsf
# based do whether LDA score and p-value meet critera.

#Read in cvs line by line, split by tab, asv seq should be index 0, P-value
#index 1 and LDA index 4.

def lefseSigFeatures (lines):
    counter = 0
    sigFeat =[]
    for line in lines:
        words= line.split('\t')
        if (abs(float(words[5])) >= 2.0 and abs(float(words[1])) <= .05):
            sigFeat.append(words[0])
            counter +=1 
    return(sigFeat,counter)
def rnaSeqSigFeatures (lines):
    counter = 0
    sigFeat =[]
    for line in lines:
        words = line.split('\t')
        if abs(float(words[4]) <= .05):
            sigFeat.append(words[0])
            counter += 1
    return(sigFeat,counter)
def metaGenomeSeqSigFeatures (lines):
    counter =0
    sigFeat =[]
    for line in lines:
        words = line.split('\t')
        if abs(float(words[1])<=.05):
               sigFeat.append(words[0])
               counter +=1
    
    return(sigFeat,counter)
def findSig(lefseResults,rnaSeqResults,metaGenomeSeqResults):  
    lefseAsvs = lefseResults[0]
    numSigLefse = lefseResults[1]
    rnaSeqAsvs = rnaSeqResults[0]
    numSigRnaSeq = rnaSeqResults[1]
    metaGenomeAsv = metaGenomeSeqResults[0]
    numSigMetaGenome = metaGenomeSeqResults[1]
    sharedAll =0
    sharedAllAsv =[]
    sharedLefseRnaSeq=0
    sharedLefseRnaSeqAsv=[]
    sharedlefseMetaGenome =0 
    sharedlefseMetaGenomeAsv =[]
    sharedRnaSeqMetaGenome =0
    sharedRnaSeqMetaGenomeAsv =[]
    for asv in lefseAsvs:
        if asv in rnaSeqAsvs:
            if asv in metaGenomeAsv:
                sharedAll +=1
                sharedAllAsv.append(asv)
            else :
                sharedLefseRnaSeq+=1
                sharedLefseRnaSeqAsv.append(asv)
        if asv in metaGenomeAsv ==True and asv in sharedAll == False:
            sharedlefseMetaGenome+= 1
            sharedlefseMetaGenomeAsv.append(asv)
    for asv in metaGenomeAsv:
        if asv in rnaSeqAsvs:
            sharedRnaSeqMetaGenome += 1
            sharedRnaSeqMetaGenomeAsv.append(asv)
        
    output =open('C:/Users/jjcol/Desktop/CeliacResearchProject/GFD healthy/Sigfeaturesall.txt','w')
    output.write(str(sharedAll) +' feature shared between RNA seq, LEFSe and metagenome seq.\n')
    for asv in sharedAllAsv:
        output.write(assignTaxonomy(asv))
    output.close()
    output =open('C:/Users/jjcol/Desktop/CeliacResearchProject/GFD healthy/SigfeaturesLefseRnaSeq.txt','w')
    output.write(str(sharedLefseRnaSeq) +' features shared only between RNA seq, LEFSe.\n')
    for asv in sharedLefseRnaSeqAsv:
        output.write(assignTaxonomy(asv))
    output.close()
    output =open('C:/Users/jjcol/Desktop/CeliacResearchProject/GFD healthy/SigFeaturesLefseMetagenomeSeq.txt','w')
    output.write(str(sharedlefseMetaGenome) +' features shared between only LEFSe and metagenome seq.\n')
    for asv in sharedlefseMetaGenomeAsv:
        output.write(assignTaxonomy(asv))
    output.close()
    output=open('C:/Users/jjcol/Desktop/CeliacResearchProject/GFD healthy/SigfeaturesRnaSeqMetaGenomeSeq.txt','w')
    output.write(str(sharedRnaSeqMetaGenome) +" Features shared between RNA seq, LEFSe and metagenome seq.\n")
    for asv in sharedRnaSeqMetaGenomeAsv:
        output.write(assignTaxonomy(asv))
def assignTaxonomy (asv):
    #change path to taxonomy file
    input= open('C:/Users/jjcol/Desktop/CeliacResearchProject/GFD healthy/taxaTableNewControls.txt')
    lines =input.readlines()
    del lines[0]
    for line in lines:
        words = line.split('\t')
        if words[0] == asv:
            return line
        
#change path to file
lefse = readIn('C:/Users/jjcol/Desktop/CeliacResearchProject/GFD healthy/lefse_de_output (3).txt')
lefseResults= lefseSigFeatures(lefse)
#change path to file
rnaSeq = readIn('C:/Users/jjcol/Desktop/CeliacResearchProject/GFD healthy/rnaseq_de.txt')
rnaSeqResults = rnaSeqSigFeatures(rnaSeq)
#change path to file 
metaGenomeSeq = readIn('C:/Users/jjcol/Desktop/CeliacResearchProject/GFD healthy/metageno_de_output.txt')
metaGenomeSeqResults = metaGenomeSeqSigFeatures(metaGenomeSeq)
findSig(lefseResults,rnaSeqResults,metaGenomeSeqResults)
