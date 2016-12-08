# -*- coding: utf-8 -*-
"""
Created on Wed May  6 12:23:02 2015

@author: legitz7
"""
#import MetricWidgets
import time, datetime
import numpy as np
import copy


#Here is the main wrapper function / API to all methods within this library, used to enable generalized thread-handling routines for
#all calls to this library made from within the GUI
def runMetric(typeName, seqList, specsList):
    #RoutineList={'CharCounter':{'func':'CharCounter', 'widget':None}   }
    if typeName not in RoutineList.keys():
        retDict={'retstatus':-1, 'message':'ERROR:  requested Metric ('+str(typeName)+') not found in master RoutineList'}
        return retDict
    else: #This is the preferred loop, where the requested Metric key was found in our RoutineList
        #Ive thought about doing this two ways, i) where we trigger production of the widget here, and
        #ii) where the widget generates the specsList that is fed straight to the function, which I think I like better
        #ii) corresponds to removal of condition below:        
        #if RoutineList[typeName]['widget'] is not None:
        #    mw=MetricWidget.MetricWidget(seqList, specsList, widget=RoutineList[typeName]['widget'] )
        retDict=RoutineList[typeName]['func'](seqList, specsList)
        retDict['seqList']=seqList    
        retDict['typeName']=typeName   
        return retDict
  
"""
OK, so charCounter is probably the most basic/fundamental Metric that we have. The idea is to produce a dictionary giving the frequency of each character combination
of a user-specified length within the *single* input sequenceObj (the function by default only operates on the first sequenceObj of a seqList passed as input)
    
There are two recent modifications to the structure of this function. 
i. For various reasons, it makes sense to structurally isolate charDicts corresponding to distinct lengths.
Therefore, we use a simpler core function ("charCounter" proper) accepting a single domainLength value, which can be called from the wrapper function accepting a range of consecutive domainLength values
We can use this to add recursion to clipping arguments fed to the outer function
 
ii. clipping arguments can be given in two distinct fashions. If a clipRange is specified, then both low and high bounds on char frequency are imposed.
Alternatively, users can specify a clipLo and/or a clipHi, each imposing the requirement that returned charDict keys occur with frequency >= clipLo, a/o < clipHi

If none of these three clipping parameters are given, no clipping will occur

###Guide to input parameters needed in the specsList:
MANDATORY
specsList['domainLength']: gives substring-length to be categorized
OPTIONAL
specsList['clipRange']:used to specify output clipping range as a two-component list 
specsList['clipLo']:used to explicitly specify lower output threshold for charDict inclusion 
specsList['clipHi']: used to explicitly specify upper output threshold for charDict inclusion
specsList['clipScale']: set this to 'pct' to specify clipBounds (via either method) as percentages rather than absolutes
UNIMPLEMENTED
specsList['forceDict']: specifies in advance a set of "valid" substrings to count
"""    
def CharCounter(seqList, specsList):
    charList=[]
    charDict={}
    timestart=time.time()
    sequence=seqList[0].Sequence
    domainLength = specsList['domainLength']
    #no need to crash the program for stupid domainLengths, filter that out and throw an Error if you see it
    if domainLength <= len(sequence):
        for i in range(len(sequence)+1-domainLength):
            if sequence[i:i+domainLength] in charList:
                charDict[sequence[i:i+domainLength]]+=1
            else:
                charList.append(sequence[i:i+domainLength])
                charDict[sequence[i:i+domainLength]]=1
        #so at this point we have a list of keys, and a dict of their values. Thus, we check the
        #clipping conditions, and filter the charDict accordingly
        charList=sorted(charList, key=lambda somechar:charDict[somechar])

        clipLo = None
        clipHi = None
        if 'clipRange' in specsList.keys():
            clipLo = specsList['clipRange'][0]
            clipHi = specsList['clipRange'][1]
        else:
            if 'clipLo' in specsList.keys():
                clipLo = specsList[clipLo]
            if 'clipHi' in specsList.keys():
                clipHi = specsList[clipHi]                       
        #We can also enter the clipping bounds as percentages of the sequence length, should we so choose
        #In that case we calculate the raw clip bounds here
            if specsList['clipScale']=='pct':
                if clipLo is not None:
                    clipLo=clipLo*(seqList[0].seqLength+1-domainLength)/100.0  
                if clipHi is not None:
                    clipHi=clipHi*(seqList[0].seqLength+1-domainLength)/100.0
        #And now we perform the filtering operation given the clipping bounds just parsed 
            if clipLo is not None:
                for echar in charDict:
                    if charDict[echar]<clipLo:
                        charList.remove(echar)
                        charDict.pop(echar, 0)
            if clipHi is not None:
                for echar in charDict:
                    if charDict[echar]>=clipHi:
                        charList.remove(echar)
                        charDict.pop(echar, 0)
        charPct={}
        for ekey in charDict.keys():
            charPct[ekey]=float(charDict[ekey])/float(seqList[0].seqLength)
        #All done, now just format and return output
         
        timestop=time.time()
        runtime=timestop-timestart
        rundate=datetime.datetime.now()    
        retDict= {'charsFound':charList, 'charFreq':charDict, 'charPct': charPct, 'typeName':'CharCounter', 'seqList':[seqList[0]], 'runtime':runtime, 'rundate':rundate}
    else:
        retDict=retDict={'retstatus':-1, 'message':'ERROR:  requested domainLength exceeds sequence length'}
    return retDict
        
"""
And here we have the outer charCounter function, which still operates on ONLY THE FIRST SEQUENCE WITHIN A seqList
but returns results for substrings of a (consecutive) range of lengths, which is useful for many more complex metrics that rely upon this method

specsList['domainRange'] needs to be a list of two bounds for character length to iterate over, i<=dR[0]<dR[1] to match "for i in range(dR[0], dR[1] )" statement


"""       
def charCounterOuter(seqList, specsList):   
    timestart=time.time()
    #The key difference in data organization between Outer and core is that all of the output dictionaries of the core function are now embedded into lists, each listitem representing a single domainLength
    charListAgg=[]
    charDictAgg=[]
    charPctAgg=[]
    #for convenience, the specsList copied to a local variable, assigned a 'domainLength' value,
    #and is passed directly to the core charCounter function for each iteration of the domainLength loop   
    
    for domainLength in range(specsList['domainRange'][0], specsList['domainRange'][1] + 1):
        localSpecs=specsList
        localSpecs['domainLength'] = domainLength
        localRetDict=CharCounter(seqList, localSpecs)
        charListAgg.append(localRetDict['charsFound'])
        charDictAgg.append(localRetDict['charsFreq'])
        charPctAgg.append(localRetDict['charPct'])
#In the future, would like to streamline the algorithm by applying forceDict based on previous loop 'charsFound' output,
#but this requires some redundant parsing of the clipping parameter inputs
        
    timestop=time.time()
    runtime=timestop-timestart
    rundate=datetime.datetime.now()  
    return {'charsFound':charListAgg, 'charFreq':charDictAgg, 'charPct': charPctAgg, 'typeName':'CharCounter', 'seqList':[seqList[0]], 'runtime':runtime, 'rundate':rundate}
            

"""
the DotPlot method is really the second "core" function of the module, since it is called by essentially all higher-level sequence comparison functions

Note for the sake of sanity that the rowSeq, which defines the equivalence-pair element specific to each row, is itself a column, and colSeq is likewise itself a row

###Guide to input parameters needed:
MANDATORY
seqList: the important thing here is that if the seqList contains 1 sequenceObj, that sequenceObj will be used for both row and column. Otherwise, seqList[0] is row and [1] col
OPTIONAL
specsList['minDisp']:used to specify minimum stretch of 'significance' (equality for default binary-equivalence weighting) nec. for recognition / storage. Defaults to 1, but obv. this will be a mess for genome-sized comparisons 
specsList['showComp']:

Long story short, default is to wmatrix mode, with binary default matrix, but 'binary' option set to False. 
specsList['binary']: when you need this to run efficiently, and are only concerned with character equality, ie no weighting matrix. This mode should be set to 'True' for most runs
specsList['wmatrix']: The alternative to binary mode. A square equivalence matrix is used, with entries given in the 'wChars' list
specsList['wchars']
specsList['dispThresh']: To use the weighting matrix, you still need to specify a "threshold"  WM value to mark a binary-yes 'dot' in the output array. Defaults to 1

UNIMPLEMENTED
"""

           
def DotPlot(seqList, specsList):
    timestart=time.time()
    #1st we assign our main variables the sequence strings rowSeq and colSeq
    rowSeqObj = seqList[0]
    rowSeq = rowSeqObj.Sequence
    if len(seqList)>1:
        colSeqObj = seqList[1]
        colSeq = colSeqObj.Sequence
    else:
        colSeqObj = rowSeqObj
        colSeq = rowSeq
    #Then we initialize a NumPy array to store the correct quantity of comparison data
    rowCount=rowSeqObj.seqLength
    colCount = colSeqObj.seqLength
    dparray=np.zeros((rowCount, colCount), dtype = int)
    #specsList optional argument parsing
    if 'dispThresh' in specsList.keys():
        dispThresh = specsList['dispThresh']
    else:
        dispThresh=1
        specsList['dispThresh']=dispThresh
        
    if 'minDisp' in specsList.keys():
        minDisp = specsList['minDisp']
    else:
        minDisp=1
        specsList['minDisp']=minDisp

    if 'showComp' in specsList.keys():
        showComp=specsList['showComp']
    else:
        showComp=False
        specsList['showComp']=showComp
# ok, an important parameter here is the weighting matrix used for scoring.
#It needs to consist of an lxl np.array, where l is the length of the alphabet. Characters are assigned to positions in the matrix via the
# wmatrixChars list (or 1-d np.array).

#All of this can be overridden by setting a 'binary' key of the specsList to True 
    #binary, accelerated option here

    binary = False
    if 'binary' in specsList.keys():
        if specsList['binary']:
            binary = True
    if binary:
        for erow in range(rowCount):
            rowVal = rowSeq[erow]
            for ecol in range(colCount):
                if colSeq[ecol] == rowVal:
                    dparray[(erow, ecol)] = 1
    #How to we want to factor the minDisp code for binary mode?
    else:
    #one will note that the weighting matrix defaults to the standard binary equivalence streamlined in the 'binary'-mode loop above
        if 'wmatrixChars' in specsList.keys():
            wmChars = specsList['wmatrixChars']
        else:
            wmChars = ['A', 'C', 'G', 'T']
        if 'wmatrix' in specsList.keys():
            wm = specsList['wmatrix']
        else:
            wm=np.diag([1]*len(wmChars))
    #filling out the raw array by performing weighting matrix lookups for each entry
        for erow in range(rowCount):
            wmrowChar = rowSeq[erow]
            wmrow = wmChars.index(wmrowChar)
            for ecol in range(colCount):
                wmcolChar = colSeq[ecol]
                wmcol = wmChars.index(wmcolChar)
                posval = wm[(wmrow, wmcol)]
                dparray[(erow, ecol)] = posval
                
        #ok now we need to handle line filtering for minDisp>1. This may be best implemented as separate helper function
        if minDisp>1:
            farray=np.zeros((rowCount, colCount))
            rarray=copy.deepcopy(farray)
        #for minDisp, there are two scans, forwards and backwards (diagonal down-right and up-left)
         #forwards scan, fills farray   
            for mrow in range(rowCount+1-minDisp):
                for mcol in range(colCount+1-minDisp):
                    if dparray[(mrow, mcol)]>=dispThresh:
                        thisl=1
                        ltermed=False
                        while not ltermed:
                            try:
                                lq=dparray[(mrow+thisl, mcol+thisl)]
                                if lq>=dispThresh:
                                    thisl+=1
                                else:
                                    ltermed=True
                            except:
                                ltermed=True
                        if thisl>=minDisp:
                            farray[(mrow, mcol)]=thisl #note that the values can often > 1. Binary mask is easy to apply tho
            
            #Backwards loop, fills rarray            
            for mrow in range(rowCount+1-minDisp):
                for mcol in range(minDisp-1, colCount):
                    if dparray[(mrow, mcol)]>=dispThresh:
                        thisl=1
                        ltermed=False
                        while not ltermed:
                            try:
                                lq=dparray[(mrow+thisl, mcol-thisl)]
                                if lq>=dispThresh:
                                    thisl+=1
                                else:
                                    ltermed=True
                            except:
                                ltermed=True
                        if thisl>=minDisp:
                            rarray[(mrow, mcol)]=thisl
        else:
        #We don't need to make these intermediate arrays if we don't have to meet significance thresolds ie MinDisp = 1
            farray=None
            rarray=None
        #and if we're showing comps too, check them as well
        #ie we are essentially repeating the routine if we need to check complementary sequences too. 
        #Should this be a separate (wmChars) call to the DotPlot routine?? 
            
        #This code may very well be redundant, comps called in by specifying approp. (actually, binary) weighting matrix
        #I wrote it again here, to simplify side-by-sides of f.,r.,c.,rc dotplots    
        if showComp:
            if minDisp>1:
                for mrow in range(rowCount+1-minDisp):
                    for mcol in range(colCount+1-minDisp):
                        if (-1*dparray[(mrow, mcol)])>=dispThresh:
                            thisl=1
                            ltermed=False
                            while not ltermed:
                                try:
                                    lq=dparray[(mrow+thisl, mcol+thisl)]
                                    if (-1*lq)>=dispThresh:
                                        thisl+=1
                                    else:
                                        ltermed=True
                                except:
                                    ltermed=True
                            if thisl>=minDisp:
                                farray[(mrow, mcol)]=-1*thisl
                for mrow in range(rowCount+1-minDisp):
                    for mcol in range(minDisp-1, colCount):
                        if (-1*dparray[(mrow, mcol)])>=dispThresh:
                            thisl=1
                            ltermed=False
                            while not ltermed:
                                try:
                                    lq=dparray[(mrow+thisl, mcol-thisl)]
                                    if (-1*lq)>=dispThresh:
                                        thisl+=1
                                    else:
                                        ltermed=True
                                except:
                                    ltermed=True
                            if thisl>=minDisp:
                                rarray[(mrow, mcol)]=-1*thisl            
                                    
    timestop=time.time()
    runtime=timestop-timestart
    rundate=datetime.datetime.now()
                
    return { 'typeName':'dotPlot', 'dparray':dparray, 'farray':farray, 'rarray':rarray, 'seqList':seqList, 'runtime':runtime, 'rundate':rundate}

    #return dparray
  
"""
OK from here on out, things get considerably more involved. Following routines include
needleman wunsch similarity
smith waterman similarity
pseudosw
summplot
modsw

Briefly review this code to see where it stands??

"""

  
    
"""
Comments on specsList: 
by default we assume that a dotplot has not already been created, and thus the function
kicks off with a dotplot of sequences to be compared.


"""    
#currently approx 350l

def needlew(seqList, specsList):
    rowSeqObj = seqList[0]
    rowSeq = rowSeqObj.Sequence
    if len(seqList)>1:
        colSeqObj = seqList[1]
        colSeq = colSeqObj.Sequence
    else:
        colSeqObj = rowSeqObj
        colSeq = rowSeq
    rowCount=rowSeqObj.seqLength
    colCount = colSeqObj.seqLength
    print rowSeq
    print colSeq
    
    hasDotPlot = False
    if 'hasDotPlot' in specsList.keys():
        if specsList['hasDotPlot']:
            hasDotPlot = True
    
    if 'wmatrixChars' in specsList.keys():
        wmChars = specsList['wmatrixChars']
    else:
        wmChars = ['A', 'C', 'G', 'T']
    if 'wmatrix' in specsList.keys():
        wm = specsList['wmatrix']
    else:
        wm=np.diag([1]*len(wmChars))    
    if hasDotPlot:
        dotPlot = specsList['dotPlot']

    else:
        dotPlot = DotPlot([rowSeqObj, colSeqObj], {'wmatrixChars':wmChars, 'wmatrix':wm})['dparray']
    print dotPlot
    if 'gapPenalty' in specsList.keys():
        gapPen=abs(specsList['gapPenalty'])
    else:
        gapPen = 2
        
    if 'returnAll' in specsList.keys():
        returnAll = specsList['returnAll']
    else:
        returnAll=True
        
    nwArray= np.zeros((rowCount + 1, colCount + 1), dtype=int)
    for i in range(1, rowCount+1):
        nwArray[(i,0)]=-gapPen*i
    for j in range(1, colCount+1):
        nwArray[(0,j)]=-gapPen*j
    #At this point, we've initialized the generator values for dynamic-programming
    #recursion, so we can just iterate through the rows to identify best m.xn. for each cell
    for m in range(1, rowCount+1):
        for n in range(1, colCount+1):
    #recall, we are always comparing three values: M(m-1,n)-g: M(m,n-1)-g: M(m-1,n-1)+W(m,n)
            val1=nwArray[(m-1,n)]-gapPen
            val2=nwArray[(m,n-1)]-gapPen
            val3=nwArray[(m-1,n-1)]+dotPlot[(m-1,n-1)]
            nwArray[(m,n)]=max(val1, val2, val3)

#ok, at this point weve made the matrix, only thing to do is to convert that to a path or paths,
#depending on whether we are in returnOne or returnAll mode...

#The sequences are always constructed in reverse. The big implementation difference is that if we want to see all
#alignments with the same scoring profile, we need to be able to spawn/fork these sequence-constructions in some way.

#to extend returnOne to the returnAll case, we have to replace the rowrSeq/colrSeq pair with
#a List of pairs, and likewise a list of nowPos positions. These three entities together define a trajectory (4, when you consider that position has 2 indices)

#so lets say we treat each fork as a 4-tuple list, which spawns adtl 4-tuple lists when degeneracies occur
#In that case, how does the logic-flow work for iterating through said list-of-lists? 
#A: it has to be a 5-tuple list, with a 'zeroed' boolean parameter for each list as well, then its just a standard while(InclusiveZeroed=False): for list in masterlist: if not zeroed: extend list

    if returnAll:
        masterList=[] #This is the list of forkable 5-tuple 'paths' ***replaced by dict
        maxscores=[] #each instance of score=maxscore gets its own forkable xxx5-tuplexxxDICT
        maxscore = 0
        for i in range(1, colCount+1):
            testval=nwArray[(rowCount, i)]
            if testval > maxscore:
                maxscores=[[rowCount, i]]
                maxscore=testval
            elif testval == maxscore:
                maxscores.append([rowCount, i])
        for j in range(rowCount-1, 0, -1):
            testval = nwArray[(j, colCount)]
            if testval > maxscore:
                maxscore = testval
                maxscores = [[j, colCount]]
            elif testval == maxscore:
                maxscores.append([j, colCount])
        #ok, now we initialize a path for every position in the maxscores list
        #masterList entry format is [rpos, cpos, rrowSeq, rcolSeq, zeroed]
        for maxpos in maxscores:
            #masterList.append([maxpos[0], maxpos[1], rowSeq[maxpos[0]-1], colSeq[maxpos[1]-1], False,maxscore])
            masterList.append({'maxRowPos': maxpos[0], 'maxColPos':maxpos[1], 'rowPos':maxpos[0], 'colPos':maxpos[1],'rRowSeq':rowSeq[maxpos[0]-1], 'rColSeq':colSeq[maxpos[1]-1], 'nowScore':maxscore, 'zeroed':False})
        allZeroed=False
        counter=0
        print masterList
        while (allZeroed == False ):
            print counter            
            counter+=1
            newList=[]
            print '---------------------------------'
            print counter
            #print masterList
            for epath in masterList:
                print epath
            for epos, eachpath in enumerate(masterList):
                newpaths=[]
                if eachpath['zeroed'] == False:
#                    valNR=nwArray[(eachpath['rowPos']-1, eachpath['colPos'] )]-gapPen
#                    valD=nwArray[(eachpath['rowPos']-1, eachpath['colPos']-1)] + dotPlot[(eachpath['rowPos']-1, eachpath['colPos']-1)]          
#                    valNC=nwArray[(eachpath['rowPos'], eachpath['colPos']-1)]-gapPen                 
                    valNR=nwArray[(eachpath['rowPos']-1, eachpath['colPos'] )]# gapPens have Already been calculated during matrix creation!
                    valD=nwArray[(eachpath['rowPos']-1, eachpath['colPos']-1)]          
                    valNC=nwArray[(eachpath['rowPos'], eachpath['colPos']-1)] 
                    valList=[valNR, valD, valNC]
                    maxval = max(valList)
                    maxList=[]
                    for vall in valList:
                        if vall == maxval:
                            maxList.append(1)
                        else:
                            maxList.append(0)
                    #maxcount = sum(maxList)

                    #now we do the for each, if= a max, make the eachpath list
                    if maxList[0]==1: #valNR
                        #newpaths.append([eachpath[0]-1, eachpath[1], eachpath[2]+rowSeq[eachpath[0]-2], eachpath[3]+'-'])
                        ep=copy.deepcopy(eachpath)                        
                        newpaths.append(ep)
                        #print 'm0'
                        #print newpaths[-1]
                        newpaths[-1]['rowPos']-=1
                        newpaths[-1]['rRowSeq']+=rowSeq[newpaths[-1]['rowPos']-1]
                        newpaths[-1]['rColSeq']+='-' 
                        newpaths[-1]['nowScore']=maxval
                        #if min(newpaths[-1][:2])<=1:
                        #print newpaths[-1]
                        if min(newpaths[-1]['rowPos'] , newpaths[-1]['colPos'])<=1:                        
                            #newpaths[-1]+=[True.maxval]
                            newpaths[-1]['zeroed'] = True
                    if maxList[1]==1: #valD
#                        newpaths.append([eachpath[0]-1, eachpath[1]-1, eachpath[2]+rowSeq[eachpath[0]-2], eachpath[3]+colSeq[eachpath[1]-2]])
#                        if min(newpaths[-1][:2])<=1:
#                            newpaths[-1]+=[True,maxval]
#                        else:
#                            newpaths[-1]+=[False,maxval]      
                        ep=copy.deepcopy(eachpath)                        
                        newpaths.append(ep)
                        #print 'm1'
                        #print newpaths[-1]
                        newpaths[-1]['rowPos']-=1
                        newpaths[-1]['colPos']-=1
                        newpaths[-1]['rRowSeq']+=rowSeq[newpaths[-1]['rowPos']-1]
                        newpaths[-1]['rColSeq']+=colSeq[newpaths[-1]['colPos']-1]
                        newpaths[-1]['nowScore']=maxval
                        #print newpaths[-1]
                        if min(newpaths[-1]['rowPos'] , newpaths[-1]['colPos'])<=1:
                            #newpaths[-1]+=[True.maxval]
                            newpaths[-1]['zeroed'] = True                   
                            
                    if maxList[2]==1: #valNC
#                        newpaths.append([eachpath[0], eachpath[1]-1, eachpath[2]+'-', eachpath[3]+colSeq[eachpath[1]-2]])
#                        if min(newpaths[-1][:2])<=1:
#                            newpaths[-1]+=[True,maxval]
#                        else:
#                            newpaths[-1]+=[False,maxval]                    
                        ep=copy.deepcopy(eachpath)                        
                        newpaths.append(ep)
                        #print 'm2'
                        #print newpaths[-1]
                        newpaths[-1]['colPos']-=1
                        newpaths[-1]['rColSeq']+=colSeq[newpaths[-1]['colPos']-1]
                        newpaths[-1]['rRowSeq']+='-' 
                        newpaths[-1]['nowScore']=maxval
                        #print newpaths[-1]
                        #if min(newpaths[-1][:2])<=1:
                        if min(newpaths[-1]['rowPos'] , newpaths[-1]['colPos'])<=1:
                            #newpaths[-1]+=[True.maxval]
                            newpaths[-1]['zeroed'] = True
                    #print 'newpaths'
                    for i in range(0, len(newpaths)):
                        #print newpaths[i]
                        newList.append(newpaths[i])
                else:
                    newList.append(eachpath) #once zeroed, a path is 'frozen' and just gets shuffled along in each loop
            #print 'newList'
           # for npt in newList:
               #print npt
            #now we need to do some culling of the newList to remove paths that have gone 'awry'. the parameter for 'elif difference<=x' is arbitrary.
            maxmax=0
            for epath in newList:
                maxmax=max(maxmax, epath['nowScore'])
            newnewList=[]
            for epath in newList:
                if epath['zeroed']:
                    newnewList.append(epath)
                elif maxmax - epath['nowScore']<=1:
                    newnewList.append(epath)
                        
            allZeroed=True
            masterList=newnewList
            #print 'masterList'
            #print masterList
            for eachpath in masterList:
                if eachpath['zeroed'] == False:
                    allZeroed=False
                    break
        #this is the completed masterList here
        #here we're adding the sequence beginnings for anything not terminating at (0,0)
        # since this is GLOBAL alignment
        print masterList
        for eachpath in masterList:
#            if eachpath[0]>1:
#                for i in range(eachpath[0]-1):
#                    eachpath[2] += rowSeq[eachpath[0]-i]
#            if eachpath[1]>0:
#                for j in range(eachpath[1]-1):
#                    eachpath[3] += colSeq[eachpath[1]-j] 
        
        #now we invert the reversed seqs to make them preety :)
        #the ldif loop should be merely prophylactic
            ldif = len(eachpath['rRowSeq'])-len(eachpath['rColSeq'])
            if ldif>0:
                eachpath['rColSeq']+='-'*ldif
            elif ldif <0:
                eachpath['rRowSeq']+='-'*abs(ldif)
                
            rowAlign=""
            colAlign=""
            for i in range(len(eachpath['rRowSeq'])):
                rowAlign += eachpath['rRowSeq'][-1-i]
            for i in range(len(eachpath['rColSeq'])):
                colAlign += eachpath['rColSeq'][-1-i]
            #now at this point we have merely the alignment, and so since this is GLOBAL we need to add
            #in any padding characters that were not included
            #first we add any tails
            rtrails = len(rowSeq)-eachpath['maxRowPos']
            if rtrails>0:
               rowAlign+=rowSeq[eachpath['maxRowPos']:] 
            ctrails = len(colSeq)-eachpath['maxColPos']
            if ctrails>0:
               colAlign+=colSeq[eachpath['maxColPos']:]
            #and now heads
            if eachpath['rowPos']>1:
                rowAlign=rowSeq[:eachpath['rowPos']]+rowAlign
                colAlign=' '*(eachpath['rowPos']-1)+colAlign
            elif eachpath['colPos']>1:
                colAlign=colSeq[:eachpath['colPos']]+colAlign
                rowAlign=' '*(eachpath['colPos']-1)+rowAlign                
            
#            rmcl=len(rowAlign)-len(colAlign)
#            if rmcl>0:
#                colAlign='-'*rmcl+colAlign
#            elif rmcl<0:
#                rowAlign='-'*abs(rmcl)+rowAlign
            compString=""
            for i in range(len(rowAlign)):
                if rowAlign[i] == colAlign[i]:
                    compString += '|'
                else:
                    compString += ' '
            #eachpath += [rowAlign, colAlign, compString]#5,6,7
            eachpath['rowFrag']=rowAlign
            eachpath['colFrag']=colAlign
            eachpath['compString']=compString
        print 'maxscore'
        print maxscore
        print '%d Alignments found: ' % (len(masterList))        
        for eachpath in masterList:
            print eachpath['rowFrag']
            print eachpath['compString']
            print eachpath['colFrag']                               
#Lets write the algorithm for the returnOne case first
    else:
        maxscore=0
        for i in range(1, colCount+1):
            testval=nwArray[(rowCount,i)]
            if testval > maxscore:
                maxscore = testval
                maxpos = [ rowCount,i]
        for j in range(rowCount, 0, -1):
            testval = nwArray[(j, colCount)]
            if testval > maxscore:
                maxscore = testval
                maxpos = [j, colCount]
        zeroed = False
        nowPos=maxpos
        rowrSeq = rowSeq[maxpos[0]-1]
        colrSeq = colSeq[maxpos[1]-1]
        while not zeroed:
            valNR=nwArray[(nowPos[0]-1, nowPos[1] )]-gapPen
            valD=nwArray[(nowPos[0]-1, nowPos[1]-1)] + dotPlot[(nowPos[0]-1, nowPos[1]-1) ]      
            valNC=nwArray[(nowPos[0], nowPos[1]-1)]-gapPen
            if max(valNR, valD, valNC) == valD:
                rowrSeq += rowSeq[nowPos[0]-2]
                colrSeq += colSeq[nowPos[1]-2]
                nowPos=[nowPos[0]-1, nowPos[1]-1]
            elif max(valNR, valD, valNC) == valNR:
                rowrSeq += rowSeq[nowPos[0]-2]
                colrSeq += '-'
                nowPos=[nowPos[0]-1, nowPos[1]]                        
            else:
                rowrSeq += '-'
                colrSeq += colSeq[nowPos[1]-2]
                nowPos=[nowPos[0], nowPos[1]-1] 
            if min(nowPos) <= 1:
                zeroed=True
        print 'nowPos'
        print nowPos
        if nowPos[0]>1:
            for i in range(nowPos[0]-1):
                rowrSeq += rowSeq[nowPos[0]-i]
        if nowPos[1]>0:
            for j in range(nowPos[1]-1):
                colrSeq += colSeq[nowPos[1]-j]                
        rowAlign=""
        colAlign=""
        for i in range(len(rowrSeq)):
            rowAlign += rowrSeq[-1-i]
        for i in range(len(colrSeq)):
            colAlign += colSeq[-1-i]        
            
        rmcl=len(rowAlign)-len(colAlign)
        if rmcl>0:
            colAlign='-'*rmcl+colAlign
        elif rmcl<0:
            rowAlign='-'*abs(rmcl)+rowAlign
        compString=""
        for i in range(len(rowAlign)):
            if rowAlign[i] == colAlign[i]:
                compString += '|'
            else:
                compString += ' '
        if maxpos[0]<(len(rowSeq)):
            rowAlign += rowSeq[maxpos[0]:]
        if maxpos[1]<(len(colSeq)):
            colAlign += colSeq[maxpos[1]:]            
        print rowAlign
        print compString
        print colAlign
        print maxscore

    return nwArray

#This implementation of Smith-Waterman is identical to the Needleman-Wunsch algorithm above, with THREE exceptions.
# 1) (most trivial/fundamental)- we now determine recursive/dp scores as max(val1, val2, val3, 0), ie things can never go below baseline in order to preserve locality
# 2) We have one more logic layer wrapping the path construction, since there is no requirement that the scores we are interested in
#reside only on the periphery of the score matrix. Thus, we require i) a boolean traversal matrix which is to be consulted in order to only initiate nonredundant path trajectories
#and ii) that the eachpath trajectory mapper iterate inward from the periphery of the matrix, shedding an outer 'layer' for each iteration
#3) a pair of new parameters threshOn and threshOff, which specify the beginning and ending conditions of a path, since we are no longer interested solely in the maximum score within the matrix

#we have also removed the 'returnOne' mode, which is practically meaningless for local alignments

#~250l

def swalpha(seqList, specsList):
    rowSeqObj = seqList[0]
    rowSeq = rowSeqObj.Sequence
    if len(seqList)>1:
        colSeqObj = seqList[1]
        colSeq = colSeqObj.Sequence
    else:
        colSeqObj = rowSeqObj
        colSeq = rowSeq
    rowCount=rowSeqObj.seqLength
    colCount = colSeqObj.seqLength
    print rowSeq
    print colSeq
    
    hasDotPlot = False
    if 'hasDotPlot' in specsList.keys():
        if specsList['hasDotPlot']:
            hasDotPlot = True
    
    if 'wmatrixChars' in specsList.keys():
        wmChars = specsList['wmatrixChars']
    else:
        wmChars = ['A', 'C', 'G', 'T']
    if 'wmatrix' in specsList.keys():
        wm = specsList['wmatrix']
    else:
        wm=np.diag([1]*len(wmChars))    
    if hasDotPlot:
        dotPlot = specsList['dotPlot']

    else:
        dotPlot = DotPlot([rowSeqObj, colSeqObj], {'wmatrixChars':wmChars, 'wmatrix':wm})['dparray']
    print dotPlot
    if 'gapPenalty' in specsList.keys():
        gapPen=abs(specsList['gapPenalty'])
    else:
        gapPen = 10

    if 'threshOn' in specsList.keys():
        threshOn = specsList['threshOn']
    else:
        threshOn = 12
        
    if 'threshOff' in specsList.keys():
        threshOff = specsList['threshOff']
    else:
        threshOff = 3            
    nwArray= np.zeros((rowCount + 1, colCount + 1), dtype=int)
    boolArray=np.ones((rowCount+1, colCount+1), dtype=bool)
    for i in range(1, rowCount+1):
        nwArray[(i,0)]=-gapPen*i
    for j in range(1, colCount+1):
        nwArray[(0,j)]=-gapPen*j
    #At this point, we've initialized the generator values for dynamic-programming
    #recursion, so we can just iterate through the rows to identify best m.xn. for each cell
    for m in range(1, rowCount+1):
        for n in range(1, colCount+1):
    #recall, we are always comparing FOUR values for Smith-Waterman: M(m-1,n)-g: M(m,n-1)-g: M(m-1,n-1)+W(m,n), AND 0
            val1=nwArray[(m-1,n)]-gapPen
            val2=nwArray[(m,n-1)]-gapPen
            val3=nwArray[(m-1,n-1)]+dotPlot[(m-1,n-1)]
            nwArray[(m,n)]=max(val1, val2, val3, 0)

#ok, at this point weve made the matrix, only thing to do is to convert that to a path or paths,

#The sequences are always constructed in reverse. The big implementation difference is that if we want to see all
#alignments with the same scoring profile, we need to be able to spawn/fork these sequence-constructions in some way.
    
    masterList=[] #This is the list of forkable 'paths' (implemented as dict for Python, struct for C)
    for z in range( rowCount, 0, -1): #this is the MASTER loop working inward from the array periphery        
        print 'z=%d'%(z)
        for i in range(colCount, 0, -1):
            testval=nwArray[(z, i)]
            if ((testval >= threshOn) and (boolArray[(rowCount,i)])) :
                boolArray[(z,i)]=False
                #threshscores.append([rowCount, i])
                thisPath={'maxRowPos': z, 'maxColPos':i, 'rowPos':z, 'colPos':i,'rRowSeq':rowSeq[z-1], 'rColSeq':colSeq[i-1], 'nowScore':testval, 'maxScore':testval,'termed':False, 'length':1}
                allTermed=False
                thisPathList=[thisPath]
                counter=0
                while (allTermed == False ):
                    counter+=1
                    print counter
                    newList=[]
                    #for epath in thisPathList:
                    #    print epath
                    maxmax=0
                    for epath in thisPathList:
                        maxmax=max(maxmax, epath['nowScore'])
                    for epos, eachpath in enumerate(thisPathList):
                        if eachpath['termed'] == False:
                            if boolArray[(eachpath['rowPos']-1, eachpath['colPos'])]:
                                valNR=nwArray[(eachpath['rowPos']-1, eachpath['colPos'] )]
                            else:
                                valNR=-1
                            if boolArray[(eachpath['rowPos']-1, eachpath['colPos']-1)]:
                                valD=nwArray[(eachpath['rowPos']-1, eachpath['colPos']-1 )]
                            else:
                                valD=-1
                            if boolArray[(eachpath['rowPos'], eachpath['colPos']-1)]:
                                valNC=nwArray[(eachpath['rowPos'], eachpath['colPos']-1 )]
                            else:
                                valNC=-1
#                            valNR=nwArray[(eachpath['rowPos']-1, eachpath['colPos'] )]   
#                            valD=nwArray[(eachpath['rowPos']-1, eachpath['colPos']-1 )]                                
#                            valNC=nwArray[(eachpath['rowPos'], eachpath['colPos']-1 )]
                            valList=[valNR, valD, valNC]
                            maxval = max(valList)
                            maxList=[]
                            for vall in valList:
                                if ((vall == maxval) and (maxval>=0)):
                                    maxList.append(1)
                                else:
                                    maxList.append(0)
                
                            newpaths=[]
                            #now we do the for each, if= a max, make the eachpath list
                            if maxList[0]==1: #valNR
                                ep=copy.deepcopy(eachpath)                        
                                newpaths.append(ep)
                                newpaths[-1]['rowPos']-=1
                                newpaths[-1]['rRowSeq']+=rowSeq[newpaths[-1]['rowPos']-1]
                                newpaths[-1]['rColSeq']+='-' 
                                newpaths[-1]['nowScore']=maxval
                                newpaths[-1]['maxScore']=max(maxval, newpaths[-1]['maxScore'])
                                newpaths[-1]['length']+=1
                                if maxval<=threshOff:                        
                                    newpaths[-1]['termed'] = True
            
                            if maxList[1]==1: #valD                
                                ep=copy.deepcopy(eachpath)                        
                                newpaths.append(ep)
                                newpaths[-1]['rowPos']-=1
                                newpaths[-1]['colPos']-=1
                                newpaths[-1]['rRowSeq']+=rowSeq[newpaths[-1]['rowPos']-1]
                                newpaths[-1]['rColSeq']+=colSeq[newpaths[-1]['colPos']-1]
                                newpaths[-1]['nowScore']=maxval
                                newpaths[-1]['maxScore']=max(maxval, newpaths[-1]['maxScore'])
                                newpaths[-1]['length']+=1
                                if maxval<=threshOff:                        
                                    newpaths[-1]['termed'] = True
                                    
                            if maxList[2]==1: #valNC                             
                                ep=copy.deepcopy(eachpath)                        
                                newpaths.append(ep)
                                newpaths[-1]['colPos']-=1
                                newpaths[-1]['rColSeq']+=colSeq[newpaths[-1]['colPos']-1]
                                newpaths[-1]['rRowSeq']+='-' 
                                newpaths[-1]['nowScore']=maxval
                                newpaths[-1]['maxScore']=max(maxval, newpaths[-1]['maxScore'])
                                newpaths[-1]['length']+=1
                                if maxval<=threshOff:                        
                                    newpaths[-1]['termed'] = True             
                            for i in range(0, len(newpaths)):
                                if maxmax-newpaths[i]['nowScore']<=1:
                                    newList.append(newpaths[i])
                                    boolArray[(newpaths[i]['rowPos'] , newpaths[i]['colPos'])]=False
                        else:
                            newList.append(eachpath) #once zeroed, a path is 'frozen' and just gets shuffled along in each loop
                    thisPathList=newList
                        #in NW, we applied a culling loop for awry paths. however, as long as everyone remains above threshold here, everyone remains valid
                        #this also seems important to avoid 'losing' paths due to cells being flagged with boolArray=False
#                        maxmax=0
#                        for epath in newList:
#                            maxmax=max(maxmax, epath['nowScore'])
#                        newnewList=[]
#                        for epath in newList:
#                            if epath['zeroed']:
#                                newnewList.append(epath)
#                            elif maxmax - epath['nowScore']<=1:
#                                newnewList.append(epath)
                                    
                    allTermed=True
                    #print 'masterList'
                    #print masterList
                    for eachpath in thisPathList:
                        if eachpath['termed'] == False:
                            allTermed=False
                            break
                for tpath in thisPathList:
                    masterList.append(tpath) #this is outside the while loop, but inside the for/if loop
                                                #thus, only completed trajectories are added here

    #this is the completed masterList here
    #here we're adding the sequence beginnings for anything not terminating at (0,0)
    # since this is GLOBAL alignment
    #print masterList
    for eachpath in masterList:

#        ldif = len(eachpath['rRowSeq'])-len(eachpath['rColSeq'])
#        if ldif>0:
#            eachpath['rColSeq']+='-'*ldif
#        elif ldif <0:
#            eachpath['rRowSeq']+='-'*abs(ldif)
            
        rowAlign=""
        colAlign=""
        for i in range(len(eachpath['rRowSeq'])):
            rowAlign += eachpath['rRowSeq'][-1-i]
        for i in range(len(eachpath['rColSeq'])):
            colAlign += eachpath['rColSeq'][-1-i]     

        compString=""
        for i in range(len(rowAlign)):
            if rowAlign[i] == colAlign[i]:
                compString += '|'
            else:
                compString += ' '
        #eachpath += [rowAlign, colAlign, compString]#5,6,7
        eachpath['rowFrag']=rowAlign
        eachpath['colFrag']=colAlign
        eachpath['compString']=compString
    if len(masterList)>0:
        #now we need to sort the list of 'print-ready' trajectory dicts
        #annoyingly, need to implement our own multi-parameter sort, since we need maxscore THEN length
        masterList=sorted(masterList, key = lambda ml: ml['maxScore'], reverse=True)
        #we subsort by refactoring masterList as a list of sublists with identical maxScores, each of which is subsorted independently and then
        #added to the dummy variable which the final masterList pointer is to be assigned to    
        nuevoList=[[masterList[0]]]
        currentmax=masterList[0]['maxScore']
        scoreCount=[[currentmax, 1]]
        for m in masterList[1:]:
            if m['maxScore']==currentmax:
                nuevoList[-1].append(m)
                scoreCount[-1][1]+=1
            else:
                currentmax=m['maxScore']
                nuevoList.append([m])
                scoreCount.append([currentmax, 1])
        dummyList=[]
    
        for spos, subList in enumerate(nuevoList):
            sortedSub=sorted(subList, key=lambda sl:sl['length'], reverse=True)
            dummyList+=sortedSub
            nuevoList[spos]=sortedSub
        masterList=dummyList
            
        for epos, escore in enumerate(scoreCount):
            print 'score=%d' % escore[0]
            print '%d Alignments found: ' % (escore[1])        
            for eachpath in nuevoList[epos]:
                print eachpath['rowFrag']
                print eachpath['compString']
                print eachpath['colFrag']               
    else:
        print 'no satisfactory alignment found'

    return nwArray


#this is my version of Smith Waterman, which uses the same dotplot matrix but more or less does
#away with indels.

#200l


def pseudosw(seqList, specsList):
    
    rowSeqObj = seqList[0]
    rowSeq = rowSeqObj.Sequence
    if len(seqList)>1:
        colSeqObj = seqList[1]
        colSeq = colSeqObj.Sequence
    else:
        colSeqObj = rowSeqObj
        colSeq = rowSeq
    rowCount=rowSeqObj.seqLength
    colCount = colSeqObj.seqLength
    print rowSeq
    print colSeq
    
    hasDotPlot = False
    if 'hasDotPlot' in specsList.keys():
        if specsList['hasDotPlot']:
            hasDotPlot = True
    
    if 'wmatrixChars' in specsList.keys():
        wmChars = specsList['wmatrixChars']
    else:
        wmChars = ['A', 'C', 'G', 'T']
    if 'wmatrix' in specsList.keys():
        wm = specsList['wmatrix']
    else:
        wm=2*np.diag([1]*len(wmChars))-np.ones((len(wmChars), len(wmChars)))    
    if hasDotPlot:
        dotPlot = specsList['dotPlot']

    else:
        dotPlot = DotPlot([rowSeqObj, colSeqObj], {'wmatrixChars':wmChars, 'wmatrix':wm})['dparray']
    print dotPlot
    if 'threshOn' in specsList.keys():
        threshOn = specsList['threshOn']
    else:
        threshOn = 8
        
    if 'threshOff' in specsList.keys():
        threshOff = specsList['threshOff']
    else:
        threshOff = 2  
# OK, so each diag of the dotplot represents an alignment. We scroll down each alignment, and whenever threshOn is met
#we look  forward and backward for threshOff, defining the boundary of our read. The alignment is then stored as a dict
#with attributes of startRowpos, startColpos, length, maxscore and sequence

#important difference is that here penalties are assessed for any mismatches, so we shouldn't run purely binary mismatches
#this is reflected in the distinct default wmatrix
    hitList=[]
    for offset in range(1-colCount, rowCount):
        if offset <=0:
            rowInit=0
            colInit=-1*offset
        else:
            rowInit=offset
            colInit=0
        alTermed=False
        nowScore=0
        rowVal=rowInit
        colVal=colInit
        threshHit=False
        while alTermed ==False:
            #print 'rowVal, colVal: (%d, %d)' % (rowVal, colVal)
            try:
                addVal=dotPlot[(rowVal, colVal)]
            except:
                alTermed=True
                break
            nowScore=max(nowScore+addVal, 0)
            if nowScore >= threshOn:
                threshHit=True
                startRowPos=rowVal
                startColPos=colVal
                startScore=nowScore #used for backward read
                
                rowRead=startRowPos
                colRead=startColPos  #dummy variables to distinguish from screening position                
                readScore=nowScore                
                maxScore=nowScore
                readLength=1
                while threshHit: #first the read forward
                    #print 'rowRead= %d, rowCount=%d' % (rowRead, rowCount)
                    if ((colRead<=(colCount-2)) and (rowRead<=rowCount-2)):
                        #print 'lets roll!'
                        rowRead+=1
                        colRead+=1
#                        try:
                        readScore+=dotPlot[(rowRead, colRead)]
#                        except:
#                            print 'error raised with rowRead=%d, colRead=%d' % (rowRead, colRead)
#                            print 'counts=(%d, %d)' % (rowCount, colCount)
#                            print 'dotplot size: (%d, %d)' % (np.shape(dotPlot)[0], np.shape(dotPlot)[1])
#                            threshHit=False
                        maxScore=max(readScore, maxScore)
                        if readScore <=threshOff:
                            rowMax=rowRead
                            threshHit=False
                        else:
                            readLength+=1                            
                    else:
                        rowMax=rowRead
                        threshHit=False
                startFound=False
                rowRead=startRowPos-1
                colRead=startColPos-1
                readScore=startScore
                while not startFound:
                    try:
                        readScore-=dotPlot[(rowRead, colRead)]
                        readLength+=1
                        if readScore <=threshOff:
                            rowFin=rowRead
                            colFin=colRead
                            startFound=True
                        else:
                            rowRead-=1
                            colRead-=1
                    except:
                        rowFin=rowRead+1
                        colFin=colRead+1
                        startFound=True
                matchDict={'rowStart':rowFin+1, 'colStart':colFin+1, 'readLength':readLength, 'maxScore':maxScore, 'offset':offset}
                hitList.append(matchDict)
                rowVal=rowMax
                colVal=rowVal+offset
            
            rowVal+=1
            colVal+=1
            if ((colVal>=(colCount-2)) or (rowVal>=rowCount-2)):
                alTermed=True
    if len(hitList)>0:
   
        #ok, now lets go ahead and grab the sequences from our hits, make a comparison String, and print
        for ehit in hitList:
            rowString=rowSeq[ehit['rowStart']:ehit['rowStart']+ehit['readLength']]
            colString=colSeq[ehit['colStart']:ehit['colStart']+ehit['readLength']]
            compString=""
            for compChar in range(ehit['readLength']):
                try:
                    rowChar=rowString[compChar]
                    colChar=colString[compChar]
                    if rowChar == colChar:
                        compString+='|'
                    else:
                        compString+=' '
                except:
                    ehit['readLength']=compChar
                    break
            lastmatch=compString.rfind('|')
            rowString=rowString[:lastmatch+1]
            colString=colString[:lastmatch+1]
            compString=compString[:lastmatch+1]
            ehit['rowString']=rowString
            ehit['colString']=colString
            ehit['compString']=compString
            ehit['readLength']=len(compString)
        
        #complex sorting loop, first by score then by length
        hitList=sorted(hitList, key=lambda d: d['maxScore'], reverse=True)
        nuevoList=[[hitList[0]]]
        currentmax=hitList[0]['maxScore']
        scoreCount=[[currentmax, 1]]
        for m in hitList[1:]:
            if m['maxScore']==currentmax:
                nuevoList[-1].append(m)
                scoreCount[-1][1]+=1
            else:
                currentmax=m['maxScore']
                nuevoList.append([m])
                scoreCount.append([currentmax, 1])
        dummyList=[]
    
        for spos, subList in enumerate(nuevoList):
            sortedSub=sorted(subList, key=lambda sl:sl['readLength'], reverse=True)
            dummyList+=sortedSub
            nuevoList[spos]=sortedSub
        hitList=dummyList 
        for epos, escore in enumerate(scoreCount):
            print 'score=%d' % escore[0]
            print '%d Alignments found: ' % (escore[1])        
            for eachpath in nuevoList[epos]:
                print 'rowPos, offset: %d, %d' %(eachpath['rowStart'], eachpath['offset'])
                print eachpath['rowString']
                print eachpath['compString']
                print eachpath['colString']         
        
    else:
        print 'no alignments found matching threshold'
        
    return hitList
        


#This (summPlot) is a really fun little tool :)
#In what regard is it different from CharCounter?
#It probable could have been implmented by calling CharCounter, but for various reasons it turned out to 
#be convenient to run a slightly different counting algorithm here
"""
SummPlot:
This produces the raw data for summPlot graphs, in the form of 
freqList: a list of substring-frequency dictionaries covering a user-specified length range,
and seqList: a list of frequency-sorted substring lists for each length.
 

#This is designed to operate on a single sequenceObj, and thus all seqList elements will be ignored except the first

###Guide to input parameters needed:
MANDATORY
seqList: the important thing here is that if the seqList contains 1 sequenceObj, that sequenceObj will be used for both row and column. Otherwise, seqList[0] is row and [1] col

OPTIONAL
specsList['minCount']:used to specify minimum threshold of substring frequency desired in output. Defaults to 2
specsList['maxLength']: this is the longest substring in which we are interested. Defaults to the lesser of (10, seqLength)
specsList['isCircular']: for analysis of plasmids, viruses, and other circular sequences. Defaults to False

specsList['genGraphics']: boolean signalling whether to run the mplGraphics.createFigSummPlot() construction routine, which edits the 'graphics' [dictionary] of specsOutputDict, 
    ultimately stored as the .graphics property of the so.Metric object
specsList['gMode']: indicates what to do with the graphic+metadata returned by createFigSummPlot()
    options include 'Live'(current default, generates live interactive pyplot sesh), 'Embedded'(PyQt subwidget), and 'Headless' (autosave gFile for later display)
specsList['gFile']: desired gFile address for headless autosave, or manual save from interactive display
specsList['pctMode']: 
UNIMPLEMENTED
                    
"""                 
def summPlot(seqList, specsList): 
    seq=seqList[0].Sequence
    seqLength=seqList[0].seqLength
    
    if 'minCount' in specsList.keys():
        minCount=specsList['minCount']
    else:
        minCount=1
        
    if 'maxLength' in specsList.keys():
        maxLength=specsList['maxLength']
    else:
        maxLength=min(10, seqLength)
        

    freqList=[]
    sortList=[]


        #generate a list of dictionaries and lists, just like charCounter    
    for q in range(maxLength):
        freqList.append({})
        sortList.append([])

    if seqLength > maxLength: #this should be like 99% of the time, unless youre being dumb
        #It turns out that if we are recording from a cyclic group the python 
        #negative-number indexing conventions become incredibly useful, and apt
        if ( (('isCircular' in specsList.keys() ) and (specsList['isCircular'] is True)) or (seqList[0].isCircular == True)):        
            for k in range(0, seqLength):
                for m in range(maxLength):
                    mchar=seq[k-m:k+1]
                    if mchar in freqList[m].keys():
                        freqList[m][mchar]+=1
                    else:
                        freqList[m][mchar]=1        
                        
        #However, for a linear sequence, the loop implementation is slightly more complicated
        #there are two consecutive for-loops that append the freqLists,
        #then one more forloop for parsing/formatting
    
        #There is a good reason why there are two distinct loops    
        #If you want a plot showing charFrequency up to maxLength, than for
        #the first (seqLength-maxLength) bases, we will be collecting substrings of all requested lengths,
        #However, as we near the end of the sequence (assuming noncircular!), then we are only able to form substrings of increasingly shorter maximum lengths,
        #and the two cases must be handled quite differently

        #The way we structure the loop is to identify, for each element of a sequence, all substrings terminating with that element, that are of acceptable length

        #so in the first for loop, we are "building up steam", starting by recording the length=1 characters,
        #and then looping over, depending on position, the correct quantity of substrings to record terminating at that position
        else:
            for i in range(maxLength):
                #record the 1-character substring
                onechar=seq[i]
                if onechar in freqList[0].keys():
                    freqList[0][onechar]+=1
                else:
                    freqList[0][onechar]=1
                #print 'freqList'
                #print freqList
                    
                #record the correct number of substrings terminating at position=i    
                for j in range(1,i+1): #(hint: it's this many...)
                    nchar=seq[i-j:i+1]
                    if nchar in freqList[j].keys():
                        freqList[j][nchar]+=1
                    else:
                        freqList[j][nchar]=1
                    #print 'freqListt'
                    #print freqList
    
    #This is the main case, where we are recording strings of all lengths <= maxLength                  
            for k in range(maxLength, seqLength):
                for m in range(maxLength):
                    mchar=seq[k-m:k+1]
                    if mchar in freqList[m].keys():
                        freqList[m][mchar]+=1
                    else:
                        freqList[m][mchar]=1
            #print 'freqList: '
            #print freqList

    #recall that the structure of the freqList we have been building is:
    # [{'A':14,'C':15,'G':11,T:12},{'AA':3,AC':2,...},{...},...]

                    
        #the sortList (which has not yet been sorted) is created after applying cutoff threshold for occurrence frequency
        for elength, edict in enumerate(freqList):
            for ekey in edict.keys():
                if edict[ekey]>=minCount:
                    sortList[elength].append([ekey, edict[ekey]])
        
        #print 'unsorted sortList'
        #print sortList    
        
        #and then each substring(len L) sublist is sorted by frequency, ie the second term of the substring entry
        #eg aiming for sortList[0] = [ ['G',15], ['A',12], ['C',11], ['T',8] ]
        

        if 'pctMode' in specsList.keys():
            pctMode=specsList['pctMode']
        else:
            pctMode=False
        
        for spos, slength in enumerate(sortList):
            

            #print 'slength'
            #print slength
              
            slength=sorted(slength, key=lambda val: val[1], reverse=True)

                #in this case we have to break out into sublists for sorting, based on previous character position
            #print 'postsort'
            #print slength

            if pctMode:
                for epos, entry in enumerate(slength):
                    entry[1]=entry[1]/float(seqLength)
                    slength[epos]=entry
                    
            sortList[spos]=slength                   

            #print 'confirm'
            #print sortList[spos]
        if 'genGraphics' in specsList.keys():
            if (specsList['genGraphics'] == False):
                genGraphics = False
            else:
                genGraphics = True
        else:
            genGraphics = True
            
        if genGraphics:
            gParamDict={}
            if 'gMode' in specsList.keys():
                gParamDict['gMode'] = specsList['gMode']
            else:
                gParamDict['gMode'] = 'Live'
            if 'gFile' in specsList.keys():
                gParamDict['gFile'] = specsList['gFile']
            if 'order' in specsList.keys():
                gParamDict['order']=specsList['order']
            else: 
                gParamDict['order']='freq'

            gParamDict['pctMode']=pctMode  
            if 'transpose' in specsList.keys():
                gParamDict['transpose']=specsList['transpose']
            else:
                gParamDict['transpose']=False
            
            import mplGraphics
            gDict=mplGraphics.createFigSummPlot(sortList, seqList, gParamDict)
        else:
            gDict = {}
    #this should be returned with freqList, since either one or both may be needed depending on how we want to display the data    
    #Recall, only sortList is a list of lists, while freqList is a list of dictionaries...
    return {'sortList':sortList, 'freqList':freqList, 'graphics':gDict}
            
        



#~250l

def modsw(seqList, specsList):

    if 'wmatrixChars' in specsList.keys():
        wmChars = specsList['wmatrixChars']
    else:
        wmChars = ['A', 'C', 'G', 'T']
    if 'wmatrix' in specsList.keys():
        wmatrix = specsList['wmatrix']
    else:
        wmatrix=np.diag([2]*len(wmChars))- np.ones((len(wmChars), len(wmChars)))
    if 'threshold' in specsList.keys():
        thresh = specsList['threshold']
    else:
        thresh = 5

    rowSeqObj = seqList[0]
    rowSeq = rowSeqObj.Sequence
    if len(seqList)>1:
        colSeqObj = seqList[1]
        colSeq = colSeqObj.Sequence
    else:
        colSeqObj = rowSeqObj
        colSeq = rowSeq
    rowCount=rowSeqObj.seqLength
    colCount = colSeqObj.seqLength
    
    mswFarray = np.zeros((rowCount, colCount))
    mswRarray = np.zeros((rowCount, colCount))
    fmatchList = []
    rmatchList = []
    allmatchList = []
    for i in range(1-colCount, 1): #offset of colSeq relative to rowSeq
        #run it fwd, then backward        
        fscore = 0
        rscore = 0
        numevals = min(rowCount, i+colCount)
        threshOn = False
        jval = -1
        jmax = -1
        kval = -1
        kmax = -1
        for j in range(numevals): #fwd first
            wmRowval = wmChars.index(rowSeq[j])
            wmColval = wmChars.index(colSeq[j-i])
            thisval = wmatrix[(wmRowval, wmColval)]
            fscore = max(0, fscore+thisval)
            mswFarray[(j, j-i)] = fscore
            if fscore >= thresh:
                if not threshOn:
                    threshOn = True
                    jval = j
                    jmax = fscore
                else:
                    jmax = max(jmax, fscore)
            else:
                if threshOn:
                    fmatchList.append([i,jval, j, jmax])
                    allmatchList.append([i,jval, j, jmax])
                    jval = -1
                    jmax = -1
                    threshOn = False
        if threshOn:
            fmatchList.append([i,jval, j, jmax])
            allmatchList.append([i,jval, j, jmax])
            jval = -1
            jmax = -1
            threshOn = False

        for k in range(numevals): #then reverse
            posnum = numevals-k-1
            wmRowval = wmChars.index(rowSeq[posnum])
            wmColval = wmChars.index(colSeq[posnum-i])
            thisval = wmatrix[(wmRowval, wmColval)]         
            rscore = max(0,rscore+thisval)
            mswRarray[(posnum, posnum-i)] = rscore
            if rscore >= thresh:
                if not threshOn:
                    threshOn = True
                    kval = posnum
                    kmax = rscore
                else:
                    kmax = max(kmax, rscore)
            else:
                if threshOn:
                    rmatchList.append([i,kval, posnum, kmax])
                    allmatchList.append([i,posnum, kval, kmax])
                    kval = -1
                    kmax = -1
                    threshOn = False
        if threshOn:
            rmatchList.append([i,kval, posnum, kmax])
            allmatchList.append([i,posnum, kval, kmax])
            kval = -1
            kmax = -1
            threshOn = False
            
    for i in range(1, rowCount):
        fscore = 0
        rscore = 0
        numevals = min(colCount, rowCount-i)
        threshOn = False
        
        for j in range(numevals): #fwd first
            wmRowval = wmChars.index(rowSeq[j-i])
            wmColval = wmChars.index(colSeq[j])
            thisval = wmatrix[(wmRowval, wmColval)]
            fscore = max(0, fscore+thisval)
            mswFarray[(j-i, j)] = fscore
            if fscore >= thresh:
                if not threshOn:
                    threshOn = True
                    jval = j
                    jmax = fscore
                else:
                    jmax = max(jmax, fscore)
            else:
                if threshOn:
                    fmatchList.append([i,jval, j, jmax])
                    allmatchList.append([i,jval, j, jmax])
                    jval = -1
                    jmax = -1
                    threshOn = False
        if threshOn:
            fmatchList.append([i,jval, j, jmax])
            allmatchList.append([i,jval, j, jmax])
            jval = -1
            jmax = -1
            threshOn = False
            
        for k in range(numevals): #then reverse
            posnum = numevals-k-1
            wmRowval = wmChars.index(rowSeq[posnum-i])
            wmColval = wmChars.index(colSeq[posnum])
            thisval = wmatrix[(wmRowval, wmColval)]         
            rscore = max(0,rscore+thisval)
            mswRarray[(posnum-i, posnum)] = rscore
            if rscore >= thresh:
                if not threshOn:
                    threshOn = True
                    kval = posnum
                    kmax = rscore
                else:
                    kmax = max(kmax, rscore)
            else:
                if threshOn:
                    rmatchList.append([i,kval, posnum, kmax])
                    allmatchList.append([i,posnum, kval, kmax])
                    kval = -1
                    kmax = -1
                    threshOn = False
        if threshOn:
            rmatchList.append([i,kval, posnum, kmax])
            allmatchList.append([i,posnum, kval, kmax])
            kval = -1
            kmax = -1
            threshOn = False
    print 'offset, jkstart, jkstop, maxscore'
#    for fm in fmatchList:
#        print fm
#    print '-----------'
#    for rm in rmatchList:
#        print rm
    sweepsize=thresh    
    allmatchList = fmatchList 
    for fm in fmatchList:
        fm[1] = max(0, fm[1]-sweepsize)
        print fm
    print '-----------'
    for rm in rmatchList:
        rm[1] += sweepsize   #watchout, this could absltly lead to OutOfBounds errors in some cases     
        allmatchList.append([rm[0], rm[2], rm[1], rm[3]])
        print rm    
    print 'allMatchList'
    allmatchList = fmatchList 
    allmatchList = sorted(allmatchList, key = lambda aml: aml[0])
    for am in allmatchList:
        print am
    contigs = []
    matchrows = []
    for match in allmatchList:
        if match[0] not in matchrows:
            matchrows.append(match[0])
    matchcounter = 0
    for matchrow in matchrows:
        thisrowMatches = []
        rightrow = True
        while rightrow:
            try:
                if allmatchList[matchcounter][0] == matchrow:
                    thisrowMatches.append(allmatchList[matchcounter])
                    matchcounter += 1
                else:
                    matchcounter += 1
                    rightrow = False
            except:
                break
        thisrowMatches = sorted(thisrowMatches, key = lambda ml: ml[1])

        for mpos, thesematches in enumerate(thisrowMatches):
            print 'mpos, thesematches'
            print mpos
            print thesematches
        
            hibound = thesematches[2]
            if mpos != len(thisrowMatches)-1:
                olobound = thisrowMatches[mpos+1][1]
                if olobound < hibound:
                    thisrowMatches[mpos+1] = [matchrow, min(thisrowMatches[mpos+1][1],thesematches[1]), max(thisrowMatches[mpos+1][2] , thesematches[2]), max(thesematches[3], thisrowMatches[mpos+1][3])]
                else:
                    contigs.append(thesematches)
            else:
                contigs.append(thesematches)

            
    
    np.set_printoptions(linewidth=200)        
#    print 'fwd\n'
#    print mswFarray
#    print 'rev\n'
#    print mswRarray
    print 'contigs'
    for contig in contigs:
        print contig   
        conlen = contig[2] - contig[1] + 1
        if contig[0] <=0:
            seq1 = colSeq[(contig[1]-contig[0]) : min(len(colSeq), contig[1]-contig[0]+conlen )]
            seq2 = rowSeq[contig[1]: min(len(rowSeq), contig[1]+conlen)]
        else:
            seq1=colSeq[contig[1]: min(len(colSeq), contig[1]+conlen)]
            seq2=rowSeq[(contig[1]+contig[0]) : min(len(rowSeq), contig[1]+contig[0]+conlen )]
        print seq1
        print seq2
        print '-------------'
    
RoutineList={'CharCounter':{'func': CharCounter, 'widget':None} , 'DotPlot':{'func': DotPlot, 'widget':None} , 'Needleman-Wunsch':{'func': needlew, 'widget':None} , 'Smith-Waterman':{'func': swalpha, 'widget':None} , 'Gapless Smith-Waterman':{'func': pseudosw, 'widget':None} , 'SummaryPlot':{'func': summPlot, 'widget':None}}    
    
    
    
    
  
    
    
    