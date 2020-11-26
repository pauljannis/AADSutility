"""
Tool for converting raw read data and detecting droplets from raw AADS trace.
Modes for peak value determination:
    -Midpoint: Middle value of event
    -Midpeak: Middle peak of event
    -Absolute: Highest point of event
    -Gradient: Middle peak determined by gradient = 0
    -Rewrite: Rewrite event trace to file to reduce filesize. No detection.

Author: Paul Zurek, pjz26@cam.ac.uk
v1.0: Combined raw detect (9c) and grad detect (1c) into one, added argparse and reduced rewrite
v1.1: Detection now in V (not 10-V). General information now also written to file
22/10/2018
Tested with python 3.6
"""

import numpy as np
from matplotlib import pyplot as plt
from random import random
import argparse

parser = argparse.ArgumentParser(description="""AADS droplet detection.                                 
                                 Author: Paul Zurek (pjz26@cam.ac.uk).
                                 Version 1.1 (22/10/2018)""")
parser.add_argument('method', help='Set detection method', default='midpoint',
                    choices=('midpoint', 'midpeak', 'absolute', 'gradient', 'rewrite'))
parser.add_argument('filename', help="input filename")
parser.add_argument('-o', '--output', help='output filename')
parser.add_argument('--verbose', action='store_true')
parser.add_argument('--manualgate', nargs=2, type=float,
                    help='set fusion gate manually (ms)')
parser.add_argument('--nogate', action='store_true',
                    help='no size gating')
parser.add_argument('--select', nargs=2, type=int,
                    help='select subset in ms to be analysed')
parser.add_argument('-v', '--version', action='version', version='1.1')

args = parser.parse_args()
v = args.verbose
out = args.output

#Some hardcoded values that you might need to change
min_event_size = 0.2  #minimal event size that is detected (ms)
max_event_size = 4    #maximal event size that is detected (ms)




#FUNCTIONS
#Load data, convert and rewrite
def loaddat(filename):
    dat = np.loadtxt(filename)
    xn, yn = zip(*dat)
    xn = np.array(xn)
    yn = np.array(yn)
    return xn, yn

#Select range of data in ms to be analysed
def select(left, right): 
    for i in range(len(xn)):
        if xn[i] < left:
            l = i
        if xn[i] < right:
            r = i
    xs = xn[l:r]
    ys = yn[l:r]
    return xs, ys

#Detect background level and set variation threshold
def background(bgv):
    av = sorted(yn)
    bg = av[int(len(av)/2)]
    bgH = bg + bgv
    bgL = bg - bgv
    return bg, bgH, bgL

#Start with and end with a nr of points at bg level (no broken droplets)
def bgstart(lst, nr):
    #Left side
    l = 0
    count = 0
    while count < nr:
        if bgL < lst[l] < bgH:
            count += 1
        else:
            count = 0
        l += 1
    #Right side
    r = len(lst) - 1
    count = 0
    while count < nr:
        if bgL < lst[r] < bgH:
            count += 1
        else:
            count = 0
        r -= 1
    ###
    return l-nr, r+nr

#Extracts events from raw data by checking if signal rises from background level and comes back to bg level
#bgnum = number of datapoints at bg level necessary for event (*f)
def extract(bgnum, f):
    frac = bgnum * f
    eventsInd = []
    up_ind = 0
    down_ind = 0
    for i in range(len(yn - bgnum)):
        if yn[i-1] > bgL > yn[i]:   #Dropping low: event start
            c = 0
            for j in range(1, bgnum):  #Check if bg level before
                if bgL < yn[i-j] < bgH:
                    c += 1
            if c > frac:
                up_ind = i #- 10
        if yn[i] < bgL < yn[i+1]:  #Going up: event complete
            c = 0
            for j in range(1, bgnum):
                if bgL < yn[i+j] < bgH:  #Check if bg level after
                    c += 1
            if c > frac:
                down_ind = i #+ 10
                size = xn[down_ind] - xn[up_ind]
                if min_event_size < size < max_event_size:  #check size appropriate, everything really should be between 0.2 and 4 ms
                    eventsInd.append([up_ind, down_ind])
    return eventsInd   #Start and end indices returned for each event

#Removes sequential indices from list
def remadjacent(lst):
    i = 0
    while i < len(lst)-1:
        if lst[i] == (lst[i+1] - 1):
            del lst[i]
        else:
            i += 1
    return lst        

#Identify peaks via null in gradient
def gradient_getpeaks(nr): 
    Gtrav = []
    g = events[nr][2]   #All gradient values
    for i in range(1, len(g)-1):
        if (g[i-1] >= 0 and g[i] <= 0) or (g[i-1] <= 0 and g[i] >= 0):
            Gtrav.append(i)
    Gtrav = remadjacent(Gtrav)
    return Gtrav

#Identify peaks via trace jiggle
def midpeak_getpeaks(nr):
    eventpeaks = []
    ys = events[nr][1]   #All y values of event
    for i in range(2, len(ys)-2):
        if ((ys[i] < ys[i-1] and ys[i] < ys[i+1]) or 
        (ys[i] < ys[i-2] and ys[i] < ys[i+2])):
            eventpeaks.append(i)
    return eventpeaks

#Plot single peak
def plotpeak(evinfo):
    nr, x, y = evinfo
    l, r = eventsIND[nr]
    xev = xn[l-10:r+10]
    yev = yn[l-10:r+10] - bg
    gev = nfg[l-10:r+10] * 5
    plt.figure(figsize=(5,10))
    plt.plot(xev, yev, "--ko")
    plt.plot(xev, gev, "--bo")
    plt.plot(x, (y - bg), "rx", markersize=15)
    plt.title("Event #%d" % nr)
    plt.ylabel("Adjusted signal [a.u.]")
    plt.xlabel("Time [ms]")
        
#Shape information of droplets for quality control in grad detect
def plotgradshapes(lst):
    shapes = [len(t) for t in lst]
    plt.figure(figsize=(10,6))
    bins = np.arange(min(shapes)-0.5, max(shapes)+0.5, 1)
    plt.hist(shapes, bins=bins)
    plt.xticks(np.arange(min(shapes), max(shapes), 1))
    plt.xlabel("Event zero points")
    plt.ylabel("Frequency")
    plt.title("Droplet shape information")
    
#Plots xn vs yn and xn vs nfg!
#Also plots list of x,y values (peak data in lst1 and lst2)
def plotfull(peakinfo, removed=[], base=10, reduced=True):
    if reduced:
        iflat1 = [range(a-base,b+base) for a,b in eventsIND]   #Convert indices into ranges
        iflat2 = [i for sub in iflat1 for i in sub]  #Flatten list of lists
        plt.figure(figsize=(10,6))
        plt.plot(xn[iflat2], yn[iflat2], "k-")
        if v: plt.plot(xn[iflat2], nfg[iflat2], "b-")
        print("%d datapoints in events (->%.1f%% for plotting)" % (len(iflat2), len(iflat2)/len(yn)*100))        
    else:
        plt.figure(figsize=(10,6))    
        plt.plot(xn, yn, "k-", lw=1)
        if v: plt.plot(xn, nfg, "b-", lw=0.5)
        print("Full dataset plotted")
    #Plot events
    nppeak = np.array(peakinfo)
    plt.plot(nppeak[:,1], nppeak[:,2], "gx")
    if len(removed) > 0:   #plot removed peaks in a different color
        rem = np.array(removed)
        plt.plot(rem[:,1], rem[:,2], "rx")
    plt.xlabel("Time [ms]")
    plt.ylabel("Detection signal [V]")
    plt.title("Raw signal")
    
#Event frequency
def freq(events):
    distances = []
    x = 0.0
    for e in events:
        if len(e[0]) > 0:
            distances.append(e[0][0] - x)
            x = e[0][0]
    distances = sorted(distances)
    med = distances[len(distances)//2]
    med = med / 1000
    hz = 1 / med
    return hz    
    
#Crop residence time by histogram density: Automatic size selection
def autogate(sizes, density):
    bins = np.arange(min(sizes)-0.0131565, max(sizes), 0.026313)   #Datapoints approx. every 26 us
    hist, bin_edges = np.histogram(sizes, bins=bins, density=True)
    #Find higest density bin
    top = 0
    for i in range(len(hist)):
        if hist[i] > hist[top]:
            top = i
    #Find left and right margin area
    bin_r = 0
    for i in range(top, len(hist)):
        if hist[i] < density:
            bin_r = i + 1
            break
    bin_l = 0
    for i in range(top, 0, -1):
        if hist[i] < density:
            bin_l = i
            break
    bin_l = bin_edges[bin_l]
    bin_r = bin_edges[bin_r]

    return bin_l, bin_r

#Write peak data to file
def output(outfile):
    out = open(outfile, "w")
    #Header information
    out.write('Dataset of %.1fs, %d events at %d Hz\n' 
              % ((xn[-1]/1000), len(eventsIND), hz))
    out.write('Final number of events extracted via %s: %d\n' 
              % (method,len(peakY)))
    out.write('Baseline voltage at:\n')
    out.write(str(bg)+'\n\nPeaks:\n')
    
    for i in range(len(residence)):
        out.write(str(residence[i]) + "\t" + str(peakY[i]) + "\n")
    out.close()
    print("droplet information written to file")
    
    
    
    
#MAIN OPERATIONS
    
filename = args.filename
#Load data
xn, yn = loaddat(filename)
print("%d datapoints loaded" % len(xn))
print("Full dataset of %.2fs" % (xn[-1]/1000))

#Select data range
sel = args.select
if sel is not None:
    xn, yn = select(sel[0], sel[1])
    print("%d datapoints selected (%d to %d ms)" % (len(xn), sel[0], sel[1]))


#Get baseline
#Increasing the value can decrease false detections for dirty samples
bg, bgH, bgL = background(0.3)  
print("Baseline is at %.2f V" % bg)


#Start at bg level (to avoid split droplets)
#nr is the number of bg level necessary for start. 40 works.
start, end = bgstart(yn, 40)
yn = yn[start:end]
xn = xn[start:end]
print("Continuous data from %d to %d" % (start, end))


#START EVENT ANALYSIS
#Extract events from dataset
eventsIND = extract(40, 0.9) #0.9 x 40 out of 40 datapoints at baselevel needed for event
print("%d events extracted" % len(eventsIND))

#Get gradient of the signal for shape control
nfg = np.gradient(yn)

#Get lists of values within start/end inidces
events = []
for e in eventsIND:
    events.append([xn[e[0]:e[1]],
                  yn[e[0]:e[1]],
                  nfg[e[0]:e[1]]])

#Calculate frequency of events
hz = freq(events)
print("Event frequency of %.1f Hz" % hz)

#DETECTION
method = args.method
if method == 'gradient':
    #Gt5: List with peak information for all 5p events:
    #[[event nr], [x of mid peak], [y of mid peak]], ...
    Gt5 = []    
    GTRAVlst = []   #List of peaks in event
    highpeak = []   #List with event information of events higher than edge
    for i in range(len(events)):
        #Get positions where the gradient is 0 (peak in trace)
        Gtrav = gradient_getpeaks(i)
        GTRAVlst.append(Gtrav)
        if len(Gtrav) == 5:  #Shape of event is as expected
            #(ID, x-value, y-value )
            Gt5.append([i, events[i][0][Gtrav[2]], events[i][1][Gtrav[2]]])
        elif 0 < len(Gtrav) < 5:
            if (bg - min(events[i][1])) > 4:  #highpeak
                highind = np.argmin(events[i][1])
                highpeak.append([i, events[i][0][highind], events[i][1][highind]])
    
    print("%d 5p shaped events" % len(Gt5))
    print("(%.1f%% of events are 5p)" % (len(Gt5)/len(eventsIND) * 100))
    #print("%d high peaks" % len(highpeak))
    #print("(%.1f%% of events are highpeak)" % (len(highpeak)/len(eventsIND) * 100))
    
    #Plot shape information for quality control
    plotgradshapes(GTRAVlst)

    #Put all events together
    eventinfo = Gt5 + highpeak

elif method == 'midpeak':
    eventinfo = []
    maxedge = 5  #highest peak instead of mid peak cutoff
    
    hmax = 0.0
    c = 0
    for i in range(len(events)):
        epeaks = midpeak_getpeaks(i)
        if len(epeaks) == 1:
            truep = epeaks[0]
        else:
            for p in epeaks:
                y = events[i][1][p]
                h = bg - y
                if h > hmax:
                    hmax = h
                    truep = p
            if hmax < maxedge:
                for p in epeaks:
                    start = events[i][0][0]
                    end = events[i][0][-1]
                    x = events[i][0][p]
                    a = (x - start) ** 2
                    b = (x - end) ** 2
                    if (a+b) < c:
                        c = a+b
                        truep = p
            c = 1000.0
            hmax = 0.0
        eventinfo.append([i, events[i][0][truep], events[i][1][truep]])
     
elif method == 'absolute':  
    eventinfo = []
    for i in range(len(events)):
        eventys = events[i][1]
        lowest = np.argmin(eventys)  #min value in event
        eventinfo.append([i, events[i][0][lowest], events[i][1][lowest]])
        
elif method == 'midpoint':  
    eventinfo = []
    for i in range(len(events)):
        eventys = events[i][1]
        mid = round(len(eventys)/2)  #middle value of event
        eventinfo.append([i, events[i][0][mid], events[i][1][mid]])

#not nicely programmed...
if method != 'rewrite':    

    #Plot some random peaks for verification
    if args.verbose:
        for i in range(0,30):
            plotpeak(eventinfo[i])
      
    #Separate eventinfo
    residence = []
    peakY = []
    for e in eventinfo:
        residence.append(events[e[0]][0][-1] - events[e[0]][0][0])
        peakY.append(e[2])
    
    #Residence histogram
    bins=np.arange(min(residence)-0.0131565, max(residence), 0.026313)
    plt.figure(figsize=(10,6))
    plt.hist(residence, bins=bins, density=True, color='black')
    plt.ylabel("Density")
    plt.xlabel("Residence time [ms]")
    plt.title("Full size histogram")
    
    
    nosizeselect = args.nogate
    manselect = args.manualgate
    removed = []
    if not nosizeselect:
        if manselect is None:
            #Size selecion: Automatic fusion removal
            l, r = autogate(residence, 0.1)   #Lower value removes less
        else:
            #Manual overwrite of crop by residence time
            l, r = manselect
        print("Detected main population at %.2f to %.2f ms"
              % (l, r))
    
        #Scatterplot
        plt.figure(figsize=(10,6))
        plt.scatter(residence, peakY, alpha=0.1)
        plt.axvline(x=l)
        plt.axvline(x=r)
        plt.ylabel("Detection signal [V]")
        plt.xlabel("Residence Time [ms]")
        plt.title("Full raw scatter with fusion gate")
        

        #Remove droplets according to size gate
        ev = len(peakY)
        for i in reversed(range(len(residence))):
            if not (l < residence[i] < r):
                removed.append(eventinfo[i])
                del residence[i]
                del peakY[i]
                del eventinfo[i]
        print("%d droplets after size selection" % len(peakY))
        print("(%.1f%% of droplets)" % (len(peakY)/ev * 100))
    
    
    #Plots trace and detected events
    reduced = not v
    plotfull(eventinfo, removed, reduced=reduced)
    
    
    #Voltage histogram
    plt.figure(figsize=(10,6))
    bins=np.arange(0,10,0.1)
    plt.hist(peakY, bins=bins)
    plt.ylabel("Frequency")
    plt.xlabel("Detection signal [V]")
    plt.title("Gated voltage histogram")
    
    
    #Gated scatterplot
    plt.figure(figsize=(10,6))
    plt.scatter(residence, peakY, alpha=0.1)
    plt.ylabel("Detection signal [V]")
    plt.xlabel("Residence Time [ms]")
    plt.title("Gated scatter plot")
    
    
    #Write peak data to output file
    if out is not None:
        output(out)
        
    plt.show()



if method == 'rewrite':
    if out is None:  #if -o is not specified the file will be overwritten!
        overwrite = input('Specify output (-o) or trace will be overwritten. Continue? (y/n): ')
        if overwrite == 'y':
            outfile = filename
        else:
            raise SystemExit('Exit: please specify filename')
    else:
        outfile = out
        
    base = 50  #Also take X values around the event for repeated detection
    iflat1 = [range(a-base,b+base) for a,b in eventsIND]   #Convert indices into ranges
    iflat2 = [i for sub in iflat1 for i in sub]  #Flatten list of lists
    f = open(outfile, "w")
    for ind in iflat2:  
        f.write(str(xn[ind]) + "\t" + str(10-yn[ind]) + "\n")
    f.close()
    print("\nTrace rewritten with eventdata:")
    print("%d datapoints in events (Data reduced to %.1f%%)" % (len(iflat2), len(iflat2)/len(yn)*100))        
 
    
    
