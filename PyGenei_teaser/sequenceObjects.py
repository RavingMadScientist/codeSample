# -*- coding: utf-8 -*-
"""
sequenceObjects Library:
Created on Tue Mar 17 18:30:58 2015


This library defines the sequenceObj, which is the core class used to represent a sequence, 
as well as various other classes that may become associated with sequenceObjs.

The library is designed to maximize the ease with which sequenceObjs can be stored and retrieved from databases.

It is important to note the distinct sequenceObj and genomeObj classes.
The two are designed to work together, with the genomeObj acting as a wrapper over the (generally quite long) genomic sequence,
allowing sequences of interest within the genome (including possibly the entire genome) to be labeled for future analysis. Unlike sequenceObjs, 
you do not run Metrics directly on a genomeObj, but rather on a desired sequenceObj obtained from the genomeObj 
     (how do we do this? is there a direct pull function?)


Current map of the library: 
we have some imports (random, sets from stdlib), and the definition of the Alphabet class, used externally for file-checking during import

Then we define two functions used in the initialization of sequenceObjs, since sequence reverses and (when applicable) complements
are stored as object properties to facilitate regex queries

 function dmap uses a dictionary to encypher a string using only characters from said dict upon initialization of sequenceObj
 function Reverse takes an input of unspecified type



Then we have a longass comment, appears to relate to original planning for the sequenceObject


Then genomeObj is defined, note brevity
Then sequenceObj, little longer

Then we define a couple companion classes, for analysis and database fun
such as the Organism class 
and a Metric (this is a nounlike container for the verbs stored in the MetricRoutines library)


Then various other obviously db-related class definitions, such as 
CodonTable, Translation, Tag

Compound, Ref, and Alphabet are defined, but essentially empty


Contrast: Sequence manipulations (functions taking as input and giving as output a single sequenceObject)
vs Metrics, which return arbitrary data, that does not typically include new sequenceObjects

To do list: identify and check MetricRoutines dependencies, is the import always necessary upon calling the library?

"""
import random
from sets import Set

# The Alphabet object is essentially just a set, generated from a list of set members
class Alphabet(object):
    def __init__(self, alpset, alpName):
        self.keys=Set(alpset)
        self.name=alpName
        
NAlphabet=Alphabet(['A','C','G','T'], 'N')
AAlphabet=Alphabet(['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'], 'A')
defaultAlps=[NAlphabet, AAlphabet]

#When sequenceObjs are initialized, we automatically store the reverse and complement (when applicable) as part of 
#the object for easier regex searching. dmap is just a convenience function to automatically replace 
# each dictionary key found within a string with the corresponding value. It generally uses the library variable compdict

compdict = {"A":"T", "C":"G", "G":"C", "T":"A" }

#friendlymode means that the function ignores all string characters not found as mapdict keys, rather than throwing an error
def dmap(inlist, mapdict, friendlymode=True):
    
    failed = False # good start
    
    try:
        mapkeys = mapdict.keys() #dict?
    except:
        print "Second argument to dmap function must be dict!"
        failed = True
        return -1
        
    if repr(type(inlist)) == "<type 'list'>":
        for pos,val in enumerate(inlist):
            if val in mapkeys:
                inlist[pos]=mapdict[val]
            else:
                if friendlymode==False:
                    print "illegal characters used in dmap; not recognized in mapdict"
                    failed = True
                    return -1
        if failed == False:
            return inlist        
            
    elif repr(type(inlist)) == "<type 'str'>":
        outstring = ""
        for val in range(0,len(inlist)):
            thischar=inlist[val]
            if thischar in mapkeys:
                outstring += mapdict[thischar]
            else:
                if friendlymode==False:
                    print "illegal characters used in dmap; not recognized in mapdict"
                    failed = True
                    return -1
                else:
                    outstring +=thischar
        if failed == False:
             return outstring
    else:
        print "first argument to dmap function must be list or string!"
        failed = True
        return -1



def Reverse(revme):
    if repr(type(revme)) == "<type 'list'>":
        try:
            melong = len(revme)
            newme=[]
            for i in range (-1, -1*(melong+1), -1):
                newme.append(revme[i])
            return newme
        except:
            print "Reverse function failed to reverse list"
            return -1
    elif repr(type(revme)) == "<type 'str'>":
        try:
            melong = len(revme)
            newme=""
            for i in range (-1, -1*(melong+1), -1):
                newme = newme + (revme[i])
            return newme
        except:
            print "Reverse function failed to reverse list"
            return -1
    else:
        print "bad argument sent to Reverse function"
        print repr(type(revme))
        print revme
        return -1
     


        
                
        
"""

This is where the SequenceObj and GenomeObj classes are stored. mandatory parameters
are all passed as arguments to the constructor function. As things stand, 
the parameters for the SequenceObj are as follows:

Mandatory:

.seqID :if object from db then this is the db PrimaryKey, if from net or 
handconstructed then this is a temp

.isGenome : is this a .visibleSequence param of a GenomeObj, or an independent entity?
.alphabet: 'N' -> nucleic, 'A' -> AA any other single-char is a custom alphabet
.Sequence : self explanatory

.source : (d=db, h=hand, n=net, f=file) (used with Dirty during database save decisions)
.
----------------------------------------------
Options + defaultAction:
.seqFile=""
.structureFiles=[]
.metrics=[] (to be implemented as fancy dict?)
.translations=[]
.parent=[]
.children=[]
.friends=[]
.tags=[]
.compounds=[]
.organism="null" (->synth is default)
.compounds=[]
.refs[]
.netsource="" (EMBL, NCBI, Uniprot, DDBJ)
.acID="" (accession ID)
---------------------------
__init()__ actions:
.Dirty calculated from .source
.rSequence : reverse of sequence
.rcSequence: if .alphabet='N', autocalculate complement sequence
.cSequence : see above
.seqLength: how long is the sequence
.masterTranslation (unsure exactly how/if this is saved to db)

Genome actions:
Mandatory:
.gFile
.vsr (view sequence range, implemented as 2-element list[], Python style)
__init__() actions:
create sequence

.seqFile: 

@author: legitz7
"""





class genomeObj(object):
    def __init__(self, gFile, vsr):
        self.gFile=gFile
        self.viewRange=vsr
        gf=open(self.gFile)
        self.seqText=gf.read()
        self.sequenceObj= ""
        gf.close()
        #companion file to raw txt of sequence, specifying feature sequence(Obj)s
        self.gFeatureFile=str(self.gFile).split('.')[0]+'_Feats'+'.txt'
        
    def addFeature(self, feature, translation=False):
        nf=open(self.gFeatureFile, 'a')
        if translation:
            nfString='TRANSLATION;'#+str(feature.range[0])+';'+str(feature.range[1])=';'+feature.tagName+'\n'
        else:
            nfString='FEATURE;'#+str(feature.range[0])+';'+str(feature.range[1])=';'+feature.tagName+'\n'
        nfString += (str(feature.range[0]) + ';')
        nfString += (str(feature.range[1]) + ';')
        nfString += (feature.tagName + '\n')
        nf.write(nfString)
        nf.close()

#annoying initialization I realize, but simplified construction enable with makeSeq function below)
#addtl optional parameters necessary for db load/store
class sequenceObj(object):
    def __init__(self,seqName,seqID, Sequence, alphabet, source, isGenome=False, myGenome="", parentID="",parentPos="",netSource="",acID="",notes="", isCircular=False):
        self.seqName=seqName
        self.seqID=seqID
        self.Sequence=Sequence.translate(None, ' ').upper() #remove whitespace before storing
        
        self.source=source
        self.isGenome=isGenome
        if self.isGenome:
            self.genome = myGenome
        self.parentID=parentID
        self.parentPos=parentPos
        self.alphabet = alphabet.upper()
        self.netSource=netSource
        self.acID=acID
        self.notes=notes
        self.isCircular = isCircular

        #Make sure illegal characters are not being used for the desired alphabet,
        #failure results in alphabet assignment of 'X'
        self.checkAlphabet()
        #ie calculate reverse and complement strings
        self.setSequence(self.Sequence)

            
        #some default conditions:::, each has a method for update
        self.metrics=[]
        self.newMetrics=[]
        self.structureFiles=[]
        self.translations=[]
        self.children=[]
        self.friends=[]
        self.tags=[]
        self.compounds=[]
        self.organism=""
        self.refs=[]
        self.seqChanged=False
        print 'Sequence initialization completed: '
        if (self.alphabet == 'A'):
            print "%s, %d aa" % (self.seqName, self.seqLength)
        else:    
            print "%s, %d bp" % (self.seqName, self.seqLength)
        ###        

            
    def checkAlphabet(self):
        legalString = True
        sl=len(self.Sequence)
        i=0        
        almap={'A':AAlphabet , 'N':NAlphabet}
        if self.alphabet in almap.keys():
            theAlphabet=almap[self.alphabet].keys
            #print theAlphabet
            for i in range(0,sl):
                if self.Sequence[i] not in theAlphabet:
                    legalString = False
        if not legalString:
            self.alphabet = 'X'

    def setSequence(self, sequence):
        self.Sequence = sequence
        sl=len(self.Sequence)
        self.seqLength = sl
        
        #for GUI preview display in tables
        if self.seqLength <= 12:
            self.preview = self.Sequence
        else:
            self.preview = self.Sequence[:6]+ '...'+self.Sequence[-6:] 
        
        #only a clean database file is not dirty
        if self.source == 'd':
            self.Dirty=False
        else:
            self.Dirty=True  

        self.rSequence = Reverse(self.Sequence)
        if self.alphabet == "N":
            self.cSequence=dmap(self.Sequence, compdict)
            self.rcSequence=Reverse(self.cSequence)
        elif self.alphabet[0]=="N":
            #this is where we would complement the extended/ambiguous code
            pass  
        else:
            self.cSequence=""
            self.rcSequence = ""
    
    #These are the class methods for editing sequence metadata    
    def addMetric(self,metric):
        self.metrics.append(metric)  
    def addstructureFile(self,sf): #raw string
        self.structureFiles.append(sf) 
    def addTranslation(self,tr):
        self.translations.append(tr) 
    def addChild(self,ch):        #sequenceObj
        self.children.append(ch) 
    def addFriend(self,friend):   # sequenceObj
        self.friends.append(friend) 
    def addTag(self,tag): 
        self.tags.append(tag) 
    def addCompound(self,compound):
        self.compounds.append(compound) 
    def addRef(self,ref):  
        self.Refs.append(ref)

#This function gives a simplified syntax for creating seqObjects from strings, 
#automatically assigning a local identifier for tracking in GUI sessions and multisequence metrics

#alp specifies whether this is a nucleic acid ('N') or amino acid ('A') sequence
def guessAlp(seq, alpList=defaultAlps):
    alpFound=False
    alpGuess = 'X'

    #search nucleotide, if found stop. else proceed to AA
    for alp in alpList:
        
        legalString = True
        sl=len(seq)
        j=0        
        for j in range(0,sl):
            if seq[j] not in alp.keys:
                #print 'illegal char for ' + alp.name + ': ' + seq[j]+'@'+str(j) 
                legalString = False
        if legalString:
            #print 'legal string for ' + alp.name
            if ((alpGuess == 'X') and (alpFound == False)):
                alpGuess = alp.name
            alpFound=True
        #print 'alpGuess = '+ alpGuess
    return alpGuess

def makeSeq(seqString, alp='N', sid='default'):
    if sid == 'default':
        qid=random.randint(0,999999999)
    else:
        qid=sid
    qs=sequenceObj('quickieSeq', qid, seqString, alp, 's')
    return qs 

def loadFromFASTA(fFile, isGenome=False):
    f=open(fFile)
    lines=f.readlines()
    f.close()
    
    header=lines[0]
    headlist=header.split('|')
    parsedDict={}
    parsedDict['seqName']=headlist[-1]
    ll=len(lines)
    seq=""
    for i in range(1,ll):
        niceline=lines[i].strip()
        seq+=niceline
    seq=seq.upper()
    parsedDict['Sequence']=seq
    parsedDict['featureTags']=[]
    seqrand=random.randint(0,999999999999)                
    seqid='n'+str(seqrand)
    alpGuess = guessAlp(seq)
    newseq=sequenceObj(parsedDict['seqName'], seqid, seq, alpGuess, source='n' )
    return newseq    

def retrieveFromDDBJ(acid, form='fasta', fileLoc='default', fileName='default', database='na'):
    import requests
    httpAddr = 'http://getentry.ddbj.nig.ac.jp/getentry'
    payload={}
    payload['database']='na'
    payload['accession_number']=acid
    payload['format']=form    
    
    r=requests.get(httpAddr,params=payload)
    if r.status_code == 200:
        print 'record retrieved: ' + acid 
#    print r.text     
    if fileName == 'default':
        fileName = acid
    if fileLoc == 'default':
        import os
        thisdir=os.getcwd()
        fileLoc=thisdir+'/data/FASTAfiles'
        if 'data' not in os.listdir('.'):
            os.mkdir(fileLoc)
        if 'FASTAfiles' not in os.listdir('./data'):
            os.mkdir(fileLoc)
    filePath=fileLoc+'/'+fileName+'.'+form
    f=open(filePath, 'w')
    f.write(r.text)
    f.close()
    
    return filePath

            
#features stored as dict {'range':[], 'name' : '' , 'featLength' : int , 'isTranslation' : bool , 'translations' : [''] }        
#def fastaParse(self,lines):
#    parsedDict={}        
#    header=lines[0]
#    headlist=header.split('|')
#    parsedDict['seqName']=headlist[-1]
#    
#    ll=len(lines)
#    seq=""
#    for i in range(1,ll):
#        niceline=lines[i].strip()
#        seq+=niceline
#    seq=str(seq.upper()) 
#    parsedDict['Sequence']=seq
#    parsedDict['featureTags']=[]
#    return parsedDict
#
#        if inMode=="File":
#            f=open(inObj)
#            lines=f.readlines()
#
#        if inStyle=='fasta':
#            defaultDict=self.fastaParse(lines)
#        self.seq=defaultDict['Sequence']
        
        
#Still in development, for the most part. Designed mostly for db load/store       
class Organism(object):
    def __init__(self,orgName,orgID,source,phyloDict={"Empire":"Synthetic"},genomeFile="",mito=False):
        self.orgName=orgName
        self.orgID=orgID
        #This basically follows the pgsql layout, {"Empire": "Kingdom":, ...}
        #quantitative analysis to be enable via upOne Dict (ie "subSpecies":"Species","Species":"Genus",...)
        #and downOne Dict, 
        self.mito=mito
        self.phyloDict = phyloDict                  
        if source=='d':
            self.Dirty=False
        else:
            self.Dirty=True

#this is a nounlike container for the verbs (functions) stored in the MetricRoutines library
#the sole MetricRoutines import in this library occurs after a Metric object is instantiated, at which time the 
# runMetric function is executed
#It makes this call in a generalized format, specifying the specific routine name, sequence(s), and run parameters
#Obviously, the metadata (specsInputDict) will be highly individual for each specific MetricRoutine 
class Metric(object):
    def __init__(self, typeName, seqList=[], specsInputDict={}, specsOutputDict={}, runName=None, seqidList=[], superdict={}, runme=True): #specsDict is params fed in, outDict is params returned

        if len(seqList)>0:        
            self.seqList=seqList
        if len(seqidList)>0:
            self.seqIDList=seqidList
        self.typeName = typeName #this is the name of the function called by the metric
        self.runName = runName #this is a user-specd title for the run
        self.specsInputDict = specsInputDict
        self.specsOutputDict = specsOutputDict      
        self.superdict=superdict
        if (specsOutputDict.keys()==[] and runme == True):        
            self.runMetric (typeName, seqList, specsInputDict)
        self.metricAddedToDB=False
        self.metricFile=""
        #self.runDate=
        #self.runtime=
    
    def runMetric(self, typeName, seqList, specsInputDict):
        import MetricRoutines
        typeList=MetricRoutines.RoutineList
        if typeName not in typeList:
            print 'ERROR: requested Metric (%s) not found in master RoutineList' % (typeName)
        else:
            self.specsOutputDict=MetricRoutines.runMetric(typeName, seqList, specsInputDict)
    
#This is primarily just an identifiable wrapper for a translation dictionary, or cdict    
class codonTable(object):
    def __init__(self, cDict={},ctName=None, ctid=None, isDefault=True):
        self.cDict = cDict
        if ctName == None:
            self.ctName='codontable'
        else:
            self.ctName = ctName
        if ctid==None:
            self.ctid='c'+str(random.randint(0,999999999999) )
        else:
            self.ctid = ctid
        self.isDefault = isDefault


#This is defined as an object for db simplicity, so that amino acid sequences remain linked to the NTs from which they derive.
#Also permits possibility of using alternative codon tables for translation, and storing results independently

#should have an alts list, for known alternate AAs, where each list element is {pos: compound}
class Translation(object):
    def __init__(self, name, parentNT, parentAA, initID="", cTable='default', tRange=None):
        self.name=name
        self.parentNT = parentNT
        self.parentAA = parentAA
        if tRange is None:
            tRange=[1,len(self.parentAA.Sequence)*3]
        self.tRange = tRange
        if len(str(initID)) < 1:
            transrand=random.randint(0,999999999999)                
            self.transID='t'+str(transrand)
        else:
            self.transID = initID
        
        if cTable == 'default':
                self.codonTable=codonTable()
                self.codonTable.cDict={"TTT":"F","TTC":"F","TTA":"L","TTG":"L","CTT":"L","CTC":"L","CTA":"L","CTG":"L",\
        "ATT":"I","ATC":"I","ATA":"I","ATG":"M","GTT":"V","GTC":"V","GTA":"V","GTG":"V","TCT":"S","TCC":"S","TCA":"S","TCG":"S",\
        "CCT":"P","CCC":"P","CCA":"P","CCG":"P","ACT":"T","ACC":"T","ACA":"T","ACG":"T","GCT":"A","GCC":"A","GCA":"A","GCG":"A",\
        "TAT":"Y","TAC":"Y","TAA":"X","TAG":"X","CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","AAT":"N","AAC":"N","AAA":"K","AAG":"K",\
        "GAT":"D","GAC":"D","GAA":"E","GAG":"E","TGT":"C","TGC":"C","TGA":"X","TGG":"W","CGT":"R","CGC":"R","CGA":"R","CGG":"R",\
        "AGT":"S","AGC":"S","AGA":"R","AGG":"R","GGT":"G","GGC":"G","GGA":"G","GGG":"G"}
        else:
            self.codonTable=codonTable
        self.parentNT.addTranslation(self)
        self.parentAA.addTranslation(self)
        
            
            
        
#designated features
class Tag(object):   
    def __init__(self, tagName, tRange='global'):
        self.tagName=tagName
        if tRange=='global':
            self.range='global'
        else:
            self.range=tRange # display-style list [1 : n] (inclusive)


#For future definition             
class Compound(object):
    def __init__(self):
        pass            

class Ref(object):  
    def __init__(self, authors, year,journal="", title=""):
        pass
        
    
        
