# -*- coding: utf-8 -*-
"""

@author: legitz7


mplGraphics is a container for the graphics objects created by MetricRoutines

It is designed to specifically encapsulate the mpl graph generation, 
returning a dictionary to be added to the so.Metric.graphics property
, and providing generalized routines for outputting these Figures:



There are generally three display modes that can be used for any metricRoutine for which "genGraphics"=True :
'Live', 'Embedded', and 'Headless'
"""

import math
import numpy as np
import matplotlib as mpl
from matplotlib.figure import Figure



#The important thing here is that we separate the code that is 
#gMode-independent (ie rectangle creation, display bounds)
#from the pyplot / PyQt / file-generation routines triggered afterwards according to gMode specified
    
#so what is conserved is that we are always going to be creating a mpl.figure.Figure, with an Axes containing the character rectangles,
#and the canvas associated with this figure is what will be used for the gMode-dependent drawing backends
    
#Unfortunately it seems like pyplot doesnt like to play nice with Figures created using base mpl routines, rather than plt.figure() for instantiation
#This is because pyplot automates the backend-dependent FigureCanvas construction, while 
#direct instantiation of mpl.figure.Figure() does not initialize any canvas by default (this is handled manually for the Qt case if gMode = 'embedded')

#Thus, it appears that in fact only the rectangle path construction is truly gMode-independent, and NOT the Figure construction

#what about custom axes instantiation??
#If we create the Axes instance explicitly, then then the graphic associated with the Axes should be reusable in all gModes
#Only problem is that the Axes object API seems to want a parent-Figure arg. to associate with right off the bat

#Perhaps it will be easiest to in fact create a direct Figure, which pyplot can scrape data from in 'Live' mode
#*** but Axes positioning within parent figure may very well be gMode-dependent, and is also a reqd instantiation argument

#Thus, the plan is:
"""
1. instantiate blank dummy Figure directly with mpl.figure.Figure()
2. instantiate Axes object directly, w req. parent set to dummy Figure
3. set Axes properties (bounds, ticks)
4. Create list of rectangle patches
5. for p in rpList: ax.add_patch(p)
-----------
gMode == 'Live': 
figOut = plt.figure()
ax.set_figure(figOut)
figOut.add_axes(ax)
ax.figure.canvas.draw()
---------
gMode == 'Embedded':
runWidget()
-------
gMode == 'Headless':
savegFile(dummyFig)

"""
def createFigSummPlot(freqList, seqList, paramDict):
    #recall typical freqList format:
    #[{'A':221, 'T':157, 'C':88,..},{..},..]
    seqLength=seqList[0].seqLength
    #Need to distinguish between freq-based v. recursive positioning
    if paramDict['order'] == 'recursive':
        posOrder='recursive'
    else:
        posOrder='freq'
    if paramDict['pctMode'] == True:
        pctMode=True
    else:
        pctMode = False
    if paramDict['transpose']== True:
        transpose=True
    else:
        transpose=False
    #defaultColorDict={'A':[0.85, 0.05, 0.0],'C':[1.0, 0.85, 0.0], 'G':[0.35, 0.35, 0.35], 'T':[0.0, 0.05, 0.50]}
    defaultColorDict={'A':[0.85, 0.05, 0.0],'C':[0.0, 0.05, 0.50], 'G':[1.0, 0.85, 0.0], 'T':[0.45, 0.45, 0.45]}
    if 'colorDict' in paramDict.keys():
        cDict=paramDict['colorDict']
    else:
        cDict=defaultColorDict
        
    gMode=paramDict['gMode']
    gDict={'gMode':gMode, 'order':posOrder, 'pctMode':pctMode, 'colorDict':cDict }
    
    
    #and actually, we should add another option here for specified sortOrders not purely alphabetical, ie to control seed positioning
    #if it's recursive, we need to resort each len(substr)=N list...
    if posOrder == 'recursive':
        #need to check for explicit ordering, 
        #else apply default (alphasorted) charOrder scraping keys from freqList[0]
        if 'charOrder' in paramDict.keys():
            charOrder=paramDict['charOrder']
        else:
            charList=[]
            for eentry in freqList[0]:
                echar=eentry[0]
                charList.append(echar)
            charOrder=sorted(charList)
        #This actually works really conveniently for proper subsorting, when you simply use alphabetical order    
        for ssLength, ssLengthList in enumerate(freqList):
            ssLengthList=sorted(ssLengthList, key = lambda entry:entry[0])
            #index-tagging each [ss,freq] entry to match formatting from freq-posOrder code below            
            for epos, entry in enumerate(ssLengthList):
                entry.append([epos])
            freqList[ssLength]=ssLengthList
            
   #In freq mode further sorting is still needed, but will be performed
    #Important to note that in following inner for loop, entry ( ['CCG',24])
    #gets appended to with a position identifier [x,y] corresponding to bot.left corner
    else: #freq mode for posOrder
        for ssLength, ssLengthList in enumerate(freqList):
            properList=[]
            if ssLength >0: #ie if not first entry
                for prevSub in freqList[ssLength-1]:
                    #print 'prevEntry'
                    #print prevSub
                    freqPos=freqList[ssLength-1].index(prevSub)
                    babyList=[]
                    for entry in ssLengthList:
                        if entry[0][:-1]==prevSub[0]:
                            entry.append([freqPos])
                            babyList.append(entry)
                    #print 'unsorted babyList'
                    #print babyList
                    babyList=sorted(babyList, key=lambda val: val[1], reverse=True )
                    #print 'complete babyList'
                    #print babyList                    
                    for entry in babyList:
                        properList.append(entry)
                freqList[ssLength]=properList
                        
                
    
#so the next thing we want to do is create the patchList,
#generate patches
    patchList=[]
    
    for ssLength, ssList in enumerate(freqList):
        #print 'ssLength: ' + str(ssLength)
        #print 'ssList:'
        #print ssList

        

        rLeft=ssLength #note, not considering rWidth
 
        prevTop=0.0
        for entry in ssList:
            #calculate upper and lower bounds
            if transpose:
                rWidth=entry[1]
                rHeight=1
                rLeft=prevTop
                rBottom=ssLength
            else:
                rWidth=1
                rHeight=entry[1]
                rLeft=ssLength
                rBottom=prevTop
            rRight=rLeft+rWidth            
            rTop=rBottom+rHeight
            patch=mpl.patches.Rectangle( (rLeft, rBottom), width=rWidth, height=rHeight)
            #apply colorMap to patches
            patchChar=entry[0][-1]            
            try:
                patchColor=cDict[patchChar]
            except:
                patchColor=[0.0,0.0,0.0] #if unknown paint it black
            patch.set_color(patchColor)
            patch.set_edgecolor([0,0,0])
            patch.set_linewidth(1)
            patchList.append(patch)
            
            if transpose:
                prevTop = rRight
            else:
                prevTop = rTop
                
            if ssLength == 0:
                #print ssLength
                #print entry
                entry.append([rLeft, rBottom])
                #print ssList
            else:
                #print ssLength
                entry[-1]=[rLeft, rBottom]
                #print entry
    gDict['patchList']=patchList

    #create dummyFig and ax, paint patches
#    dFig = Figure()
#    ax = mpl.axes.Axes(dFig, [0,0,1.0,1.0]) # left, bottom, width, height : Axes position on its dummy parent Figure
#    ax.set_xlim(left=0, right=len(freqList))   
#    ax.set_xbound(lower=0, upper=len(freqList))
#    ax.set_xticks(range(0,len(freqList)))
#    #I believe this is this the first point at which the pct v. raw_count decision matters. 
#    #but not nec. to implement immediately...
#    if paramDict['pctMode'] == True:
#        ax.set_ylim(bottom=0, top=1.0)
#        ax.set_ybound(lower=0, upper=1.0)
#        ax.set_yticks(range(0, 1.05, 0.1))  
#    else:
#        ax.set_ylim(bottom=0, top=seqList[0].seqLength)
#        ax.set_ybound(lower=0, upper=seqList[0].seqLength)
#        #we need a more complex tick scheme here..
#        slog=math.log10(seqLength)
#        naiveMarkSpacing=10**(math.floor(slog))        
#        naiveMarkCount= 1+math.floor( (seqLength/naiveMarkSpacing) )        
#        if naiveMarkCount>=5:
#            markSpacing=naiveMarkSpacing
#        elif naiveMarkCount <=2:
#            markSpacing=naiveMarkSpacing/5
#        else:
#            markSpacing=naiveMarkSpacing/2
#        ax.set_yticks(range(0, seqLength, int(markSpacing)))
#        
#    for patch in patchList:
#        ax.add_patch(patch)        
#    dFig.add_axes(ax)

    #load axes and patches to gDict
#    gDict['Axes']= ax

#and then turns out easiest way to do this is to immediately break out into
#gMode-dependent loops for figures and axes 
    
    #gMode handling
    gMode=paramDict['gMode']
    if ((gMode == 'Live') or (gMode=='Headless')): #Use same construction for Headless, just use .savefig(routine at the end)
        import matplotlib.pyplot as plt
        figOut=plt.figure()
        nax=figOut.add_subplot(111)
        for patch in patchList:
            nax.add_patch(patch)
        if transpose:
            if pctMode:
                nax.set_xlim(left=0, right=1.0)
            else:
                nax.set_xlim(left=0, right=seqLength)
            nax.set_ylim(bottom=0, top=len(freqList))  
            nax.yaxis.set_ticks(np.arange(0,len(freqList), 1))
        else:
            if pctMode:
                nax.set_ylim(bottom=0, top=1.0)            
            else:
                nax.set_ylim(bottom=0, top=seqLength)

            nax.set_xlim(left=0, right=len(freqList))
            nax.xaxis.set_ticks(np.arange(0,len(freqList), 1))
   

#        red_patch = mpl.patches.Patch(facecolor='red', linewidth=1)
#        green_patch = mpl.patches.Patch(facecolor='green', linewidth=1)
        #nax.legend(['A Legend'])
        #nax.legend(handles=[red_patch])
        cPatches=[]
        cChars=[]
        for eentry in freqList[0]:
            echar=eentry[0]
            cChars.append(echar)
            cColor=cDict[echar]
            cPatch=mpl.patches.Patch(facecolor=cColor, linewidth=1)
            cPatches.append(cPatch)
        cPatches.reverse()
        cChars.reverse()
        patchTuple=tuple(cPatches)
        charTuple=tuple(cChars)
        figOut.gca().legend(patchTuple,charTuple)
        if 'title' in paramDict.keys():
            title=paramDict['title']
        else:
            title='SummaryPlot:' + seqList[0].seqName
        plt.title(title, fontdict={'fontsize':8})
        #plt.text(0,0,'Hello')
#        legPatches=[]        
#        for echar in sorted(cDict.keys(), reverse=True):
#            patch=mpl.patches.Patch(color=cDict[echar], label=echar)
#            legPatches.append(patch)                             
#        plt.legend(handles=legPatches)
        #nax.add_artist(legend)
        #ax.set_figure(figOut)
        #figOut.add_axes(ax)
        gDict['Figure']=figOut
        gDict['Axes']=nax
        nax.figure.canvas.draw()
        if ((gMode=='Headless') and ('gFile' in paramDict.keys())):
            gFileName=paramDict['gFile']
            gDict['gFile']=gFileName            
            print 'saving fig to '+gFileName
            plt.savefig(gFileName)
            plt.close()
    elif gMode == 'Embedded':
        pass
    else:#headless, should just be a call to mplg.saveFigure()
        pass    
    return gDict
