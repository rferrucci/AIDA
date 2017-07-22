#!/usr/bin/env python

""" AIDA now accepts arlequin format for diploid data. Can use diploid-arlequin, diploid-excel, haploid, frequency
Data separation can be whitespace or none. 0 setup as own distance class. calculates confidence intervals.
Can provide own bins, or choose equal numbers or distance for bins. Adding option for microsattelites
program will dummy code microsattelites. Need to figure out if sequence option is still in the program
need to re-write for microsattlites for counting number of polymorphic sites.
"""

from numpy import resize, square, transpose, array, mean, sum,power, multiply,std, var, zeros, float32, int8
#from numpy import *

#from matplotlib.pylab import figure, plot, title, show, xlabel, ylabel, figure, savefig
from random import random
#from numpy import *     #resize, transpose
#from numpy import resize, transpose, array, mean, sum, random, zeros
from numpy.ma import masked_values

#from rrf import arp2excel, AIDAfreq, GreatCircle, Cartesian

import sys
from math import sqrt, radians, atan2, ceil
from logging import FileHandler, info
from logging import warning, basicConfig, DEBUG

####--------------------------------------------------------------------####
####---------------------default parameters-----------------------------####

# default values of input parameters:
# c = distance measure (g = great circle,c = network); i = input files (e = excel, a = arelquin)
# f = data format (h = haploid, d = diploid excel, f = population frequency, mi,ms = microsattelite IMM,SMM);
# p =  permutations; b = bins; ci = indices (b = both, g = Geary, m = Moron)
# d = distance class (d = equal distance, n = equal numbers, u = specified)
# ci = confidence limits (d = distance class; s = single-set; n = none)
# m = missing data; s = spacing (w = whitespace, n = none)

c = 'g'; i = 'e'; f = 'h'; p = 1000; d = 'n'; m = 'n'; s = 'w'; ci = 'd';
indice = 'm'

coordinates = c; format = f; distance_class = d; missing = m; space = s;
input = i; numDistClasses = 3;

permutations = p

#output= open('sylvia.out','w')

allelesfile = sys.argv[1];      del sys.argv[1]
coordinatesfile = sys.argv[1];      del sys.argv[1]
output= sys.argv[1];      del sys.argv[1]

#output = open(outfile, 'w')
out = output.split('.')
efile = '%s.error' %out[0]

efile = open(efile,'w')
sys.stderr = efile

x = allelesfile.split('.')
logfile = '%s.log' %out[0]

basicConfig(level=DEBUG,
    format='%(asctime)s %(levelname)s %(message)s',
    filename=logfile,
    filemode='w')

###-----------------options-------------------------###

from optparse import OptionParser
parser = OptionParser()

parser.add_option('-c', dest='c')
parser.add_option('-i')
parser.add_option('-f')
parser.add_option('-p',type = 'int')
parser.add_option('-d')
parser.add_option('-b')
parser.add_option('-m')
parser.add_option('-s')
parser.add_option('-a')
parser.add_option('-z')

(options,args) = parser.parse_args()


if options.c:
	coordinates = options.c
if options.i:
	input = options.i
if options.f:
	format = options.f
	if format == 'f':
		x = options['f'].split('=')
		allelefile = x[0]
if options.d:
	distance_class = options.d
if options.p:
	p = options.p
	permutations = p
if options.b: 
        bins = []
        if distance_class == 'u':
            i = 1
            x = options.b.split('=')
	    y = x[-1].split(':')
	    numDistClasses = int(x[0])
            bins = map(eval,y)
	else:
	    numDistClasses = int(options.b)

if options.z:
    	ci = options.z
if options.m:
	missing = options.m
if options.s:
	space = options.s     
if options.a:
	indice = options.a

###-------------------------------------------------###
###-------------------------------------------------###
###-------------------------------------------------###

####--------------------------------------------------------------------####
###-------------------data input----------------------------------------####

if missing == '?':
    #missing = -9
    from re import compile, sub, search
    p = compile('\?')

freqs=[];

info('importing data')

coFile = {}

pops = []

if coordinates == 'u':
    distvec = []
    for i in open(coordinatesfile):
        x = i.split()
	pops.append(x[0])
	distvec += (map(eval,x[1:]))

else:
    for i in open(coordinatesfile):
	    x = i.split()
	    pop = x[0].strip()
	    pops.append(pop)
    
	    coFile[pop] = map(eval,x[1:3])

datafile = []

for i in open(allelesfile):
    if missing == '?':
        if search('\?', i):
            i = p.sub('-9',i)
    datafile.append(i)

if input == 'e':
    data = []

    for i in datafile:     #this method purported to be quicker and easier on memory than "readlines".

        j = i.strip()
        x = j.split()
        line = x[:2]

        if space == 'n':
            line += map(eval,list(x[2]))
        elif space == 'w':
            line += map(eval,x[2:])
        #PopCount[x[1]] += 1

        data.append(line)

elif input == 'a':
    from rrf import arp2excel
    data = arp2excel(datafile, space, format)

def DummyIMM(G, alleles):
        #codes pair of alleles for homozygote presence (2) or absence (0) or hetero (1)
        #G is data dictionary, peeps is list of samples in datadict, to control order
        #of input
        dummydata = []
        #o = open('GCdummy.txt','w')
        n = 1
        for i in range(len(G)):
            b = []
            for j in range(len(G[i])):
                for a in alleles[j]:
                    if G[i][j] == -9:
                        b.append(-9)
                    elif G[i][j] == a:
                        b.append(1)
		    else:
			b.append(0)

            n += 1
            dummydata.append(b)
        #o.close()
        return dummydata

def DummySMM(G, alleles):
    dummydata = []
    for i in range(len(G)):
        s = ''
	for j in range(len(G[i])):
		r = G[i][j] - alleles[j][0] + 1
		z = alleles[j][1] - G[i][j]
		s += '1' * r + '0' * z

	dummydata.append(map(eval,list(s)))

    return dummydata

if format == 'd':
    #from rrf import indfreq
    from rrf import dip2hap
    #freqs = indfreq(data)
    data, freqs = dip2hap(data)

elif format == 'h' or format == 's':
    for line in data:
        freqs.append(line[2:])

elif format == 'f':
    from rrf import AIDAfreq
    hapfile = []
    for i in open(allelefile):
        hapfile.append(i)

    data = AIDAfreq(data,hapfile)

    del hapfile
    for i in data:
        freqs.append(i[2:])

elif format == 'mi':   #change to mi for infinite alleles
    #yr means using dummy alleles but the data is still in raw, i.e., microsattelite form
    #yd means input is already in dummy allele form
    for line in data:
        freqs.append(line[2:])

    alleles = []

    #k = 0

    A = [i for i in freqs]

    A = transpose(A)
    for i in range(len(A)):
        d = [j for j in set(A[i]) if j != -9]
        alleles += [d]

    freqs = DummyIMM(freqs, alleles)
    del alleles; del A; del d; del DummyIMM

elif format == 'ms':
    #yr means using dummy alleles but the data is still in raw, i.e., microsattelite form
    #yd means input is already in dummy allele form
    for line in data:
        freqs.append(line[2:])

    alleles = []

    #k = 0

    A = [i for i in freqs]

    A = transpose(A)
    for i in range(len(A)):
        d = [j for j in set(A[i]) if j != -9]
        a = min(d), max(d)
        alleles.append(a)

    freqs = DummySMM(freqs, alleles)
    del alleles; del A; del d; del DummySMM

else: print "ERROR: datatype not known"

#coPairs = [coFile[i[1]] for i in data]

#for i in data:
#    print i
#    pop = i[1]
#    coPairs.append(coFile[i[1]])

####-----------do some stuff with the data---------------####

if missing == -9:
    #from numpy.core.ma import masked_values
    missing = -9
    freqs = masked_values(freqs,missing)
elif missing == 'none':
    freqs = array(freqs,dtype=float32)	#convert to an array, may do this later as arrays are inflexible

npops,ngenes = len(freqs), len(freqs[1])
n =  npops

U = resize(mean(freqs,0),(npops,ngenes))		#create row vector of means
freqs = freqs - U	        #matrix of deviations from the mean
transfreqs = transpose(freqs)	#transpose of the previous matrix

del U


if format == 's':
    #from numpy import var
    transfreqs = array([x for x in transfreqs if var(x) != 0],dtype=float32)
    freqs = transpose(transfreqs)
    ngenes = len(transfreqs)
###------------------distance vector setup----------------------------###
info('calculating coordinate distances')

def Cartesian(coFile, pops):
    from math import sqrt
    distvec = [sqrt((coFile[i][0]-coFile[j][0])**2 + (coFile[i][1]-coFile[j][1])**2) for i in coFile for j in coFile]

    return distvec
    #dist = resize(distvec,(npops,npops))

def GreatCircle(coFile,pops):

    from math import sin, cos, atan2, sqrt, pi, radians

    #i = 0

    newCoors = [map(radians,coFile[p]) for p in pops]

    distvec = []
    distvec = [atan2(sqrt ( (cos(j[0]) * sin(abs((i[1]-j[1]))))**2
    	+ (cos(i[0])*sin(j[0]) - sin(i[0])*cos(j[0])*cos(abs((i[1]-j[1]))))**2),
    	sin(i[0])*sin(j[0]) + cos(i[0])*cos(j[0])*cos(i[1]-j[1]))
    	* 6367.65 for i in newCoors for j in newCoors]

    #for i in range(len(data)):
    #	for j in range(len(data)):
    #		dist[i][j] = distvec[pops.index(data[i][1])*len(pops) + pops.index(data[j][1])]

    return distvec


def degrees(d,m,s):
	s = m*60 + s
	m = s/3600.
	return d + m

if coordinates == 'c':
    #dist = zeros((len(freqs),len(freqs)),dtype=float32)

    #CoPairs = Cartesian(coFile,pops)
    distvec = Cartesian(coFile,pops)
    del Cartesian

    #distvec = Cartesian(coPairs)
    #for i in range(len(data)):
    #	for j in range(len(data)):
    #		dist[i][j] = CoPairs[pops.index(data[i][1])*len(pops) + pops.index(data[j][1])]     #[r + lp*l]

elif coordinates == 'g':
    #CoPairs, dist = GreatCircle(coFile,pops)
    distvec = GreatCircle(coFile,pops)

    del GreatCircle

dist = zeros((len(freqs),len(freqs)),dtype=float32)


for i in range(len(data)):
    for j in range(len(data)):
        #dist[i][j] = CoPairs[pops.index(data[i][1])*len(pops) + pops.index(data[j][1])]     #[r + lp*l]
        dist[i][j] = distvec[pops.index(data[i][1])*len(pops) + pops.index(data[j][1])]     #[r + lp*l]
    #distvec = []
    #i = 0

    #[r + lp*l]
from numpy import ravel


CoPairs = distvec
#elif coordinates == 'u':


#psyco.bind(GreatCircle)


    #for l in range(lp):
   #	for c in range(PopCount[pops[l]]):
#		for r in range(lp):
#			#print l,r,lp*l,CoPairs[r + lp*l]
#			#distvec += [CoPairs[r + lp*l] for i in range(PopCount[pops[r]])]
   #     		for i in range(PopCount[pops[r]]:
   #     			distvec[c][]

#npops = len(freqs)
#dist = resize(distvec,(npops,npops))
#print "dist vec finished"
#print dist

###--------------------bin time-----------------------------------###
info('set-up bins')

minvec = []

#for i in CoPairs:
for i in distvec:
        if i == 0: continue
        else: minvec.append(i)

mindist, maxdist = min(minvec), max(distvec)

upper = ceil(maxdist)

del datafile;

if distance_class == 'n':
    bins = [0]
    minvec.sort()
    classSize = len(minvec)/(numDistClasses-1)
    for i in range(0,len(minvec),classSize):
        bins.append(minvec[i])

    if len(bins) == numDistClasses:
    	bins.append(maxdist+0.01)
    elif len(bins) >= numDistClasses:
    	bins[-1] = maxdist+0.01

elif distance_class == 'u':
    bins.append(maxdist+0.01)

elif distance_class == 'd':
    "bin numeric vector"

    mi,ma=mindist, maxdist+0.01
    inc=(ma-mi)/float(numDistClasses-1)

    bins=[0,mi]
    for i in range(numDistClasses-1):
	bins.append(mi+(i+1)*inc)
	#bins.append(inc * (i+1))

    bins[-1] = ma

classDistance = []      #create list to determine average distance of distance class


for i in range(len(bins)-1):
    j = i+1
    classDistance.append((bins[i]+bins[j])/2)


###-----------------------determine distance class matrices----------###
info('determine distance class matrices')

binlength = len(bins)

DistanceClasses = []    #for list of distanceclass matrices (i.e., Wij)
distanceClass = []      #for list of distance classes with data in distanceclass



saveout = sys.stdout
sys.stdout = open(output, 'w')

def DistanceClassSetup():
	for r in range(0,binlength-1):		#set up W matrices for distance classes
	    distanceclass = r

	    #Wij = zeros((len(dist),len(dist[0])),dtype=int8)
	    #j = i+1

	    def if_else(condition, a, b):
		if condition:
		    return a
		else:
		    return b

	    n = 0
	    if distanceclass == 0:
		Wi = map(lambda x: if_else(x==0, 1, 0), CoPairs)
	    else:
		Wi = map(lambda x: if_else(x >= bins[r] and x < bins[r+1],1,0), CoPairs)

		#W.append(Wi)

	    Wij = zeros((len(freqs),len(freqs)),dtype=int8)

	    for i in range(len(freqs)):
		for j in range(len(freqs)):
			if i == j: Wij[i][j] = 0
			else:
				Wij[i][j] = Wi[pops.index(data[i][1])*len(pops) + pops.index(data[j][1])]     #[r + lp*l]

	    #Wij = array(Wij)
	    if sum(Wij) == 0:
		print "no data in distance class", distanceclass
	    else:
		distanceClass.append(i)
		DistanceClasses.append(Wij)

#psyco.bind(DistanceClassSetup)
DistanceClassSetup(); del DistanceClassSetup

del CoPairs; del data;

"""
for i in range(0,binlength-1):		#set up W matrices for distance classes
    distanceclass = i
    Wij = []

    #Wij = zeros((len(dist),len(dist[0])),dtype=int8)
    j = i+1
    s = 0
    print bins[i], bins[j]

    for popi in dist:

        Wi = []
        Wbin = []

        def if_else(condition, a, b):
            if condition:
                return a
            else:
                return b

        if distanceclass == 0:
            Wi = array(map(lambda x: if_else(x==0, 1, 0), popi),dtype = int8)
            Wi[s] = 0
        else:
            Wi = array(map(lambda x: if_else(x >= bins[i] and x < bins[j],1,0), popi),dtype = int8)
            Wi[s] = 0

        Wij.append(Wi)
        s += 1

    Wij = array(Wij)
    if sum(Wij) == 0:
        print "no data in distance class", distanceclass
    else:
        distanceClass.append(i)
        DistanceClasses.append(Wij)
"""

n =  npops

del dist

info('begin output')

###------------------initial output-------------------------###
print "Coordinates file ==========================>", coordinatesfile
print "Alleles file ==============================>", allelesfile
print "Number of polymorphic sites ===============>", ngenes
print "Number of sequences =======================>", npops
print "Minimum distance between populations ======>", '%.2f' % mindist
print "Maximum distance between populations ======>", '%.2f' % maxdist
print "for great circle distance calculations, minimum and maximum distances are in kilometers"
###---------------------------------------------------------###

freqs2 = freqs**2
freqs2 = freqs2.sum()

###------------------AIDA time-----------------------------------------###

n = npops

def AIDA(transfreqs, Wij, freqs2=freqs2,npops = npops):
    Num = (sum(Wij))*freqs2        #begin calculating numerator for Moran's II
    n = npops
    m = len(freqs[0])
    def Moran(allele):
        mat = resize(allele,(npops,npops))
        #tmat = transpose(mat)

        if missing == -9:
            mat = masked_values(mat,missing)

    	WCP = multiply(mat,transpose(mat))*Wij
        return WCP.sum()


    def Geary(allele):
    	mat = resize(allele,(npops,npops))

    	if missing == -9:
            mat = masked_values(mat,missing)

    	diff = (mat - transpose(mat))**2

        #Wdiff = diff * Wij
        #return Wdiff.sum()
        return (diff * Wij).sum()

    if indice == 'b':
        II = n*sum([Moran(x.copy()) for x in transfreqs])/Num
        cc = (n - 1)* sum([Geary(x.copy()) for x in transfreqs])/(2 * Num)
        return (II, cc)

    elif indice == 'm':
        II = n*sum([Moran(x.copy()) for x in transfreqs])/Num

        #II = n*sum([sum([freqs[i][k] * freqs[j][k] for k in range(m)])* Wij[i][j]
        #	for i in range(n) for j in range(n)])/Num

        return II

    elif indice == 'g':
        cc = (n - 1)* sum([Geary(x.copy()) for x in transfreqs])/(2 * Num)
        return cc

#if indice == 'g':
#    del Moran
#elif indice == 'm':
#    del Geary

if indice == 'b':
    print "\nCLASS LIMITS\t  II\t  cc\tW\n >=     <"
elif indice == 'm':
    print "\nCLASS LIMITS\t  II\tW\n >=     <"
elif indice == 'g':
    print "\nCLASS LIMITS\t  cc\tW\n >=     <"

MoranII = []; Gearycc = []; DistanceClass = []

i = 0

info('begin AIDA')


for Wij in DistanceClasses:
    j = distanceClass[i]

    dc = 'distance class ' + str(i)
    info(dc)

    k = j+1

    if indice == 'b':
        (II,cc) = AIDA(transfreqs,Wij)
    elif indice == 'm':
        II = AIDA(transfreqs,Wij)
    elif indice == 'g':
        cc = AIDA(transfreqs,Wij)

    #(II,cc) = AIDA(transfreqs,Wij)
    w = sum(Wij)/2

    #if II == 0:
    #       continue

    #MoranII.append(float(II))
    #Gearycc.append(float(cc))

    DistanceClass.append(distanceClass[i]+1)

    if indice == 'b':
        if bins[i] == 0:
            print '  %.2f  %.2f    %.4f  %.4f  %d' % (bins[i],bins[i],II,cc,w)
        else:
            print '%-4.2f  %-4.2f    %.4f  %.4f  %d' % (bins[i],bins[i+1],II,cc,w)
        i += 1
    elif indice == 'm':
        if bins[i] == 0:
            print '%  .2f  %  .2f    %.4f  %d' % (bins[i],bins[i],II,w)
        else:
            print '%  .2f  %  .2f    %.4f  %d' % (bins[i],bins[i+1],II,w)
        i += 1
    elif indice == 'g':
        if bins[i] == 0:
            print '%  .2f  %  .2f    %.4f  %d' % (bins[j],bins[i],cc,w)
        else:
            print '%  .2f  %  .2f    %.4f  %d' % (bins[i],bins[i+1],cc,w)

        i += 1

    #if bins[j] == 0:
    #    print '%.2f   %.2f\t%.4f  %.4f  %d' % (bins[j],bins[j],II,cc,w)
    #else:
    #    print '%.2f   %.2f\t%.4f  %.4f  %d' % (bins[j],bins[k],II,cc,w)
    #i += 1

print "\n"
b2 = (sum(freqs ** 4))/((sum(freqs**2))**2)
###----------------confidence intervals----------------------###

def scoreatpercentile(dc, perc):
    ind = int(perc * p)
    return dc[ind]

if ci == 'd':
    #from random import random
    from numpy.random import shuffle
    MoranII_ci = []; Gearycc_ci = []

    info('begin confidence intervals')

    print "PERMUTATIONS =", p, "\n"
    for dc in xrange(0,len(distanceClass)):
        j = str(dc)
        MoranII_ci.append([])
        Gearycc_ci.append([])

    for i in xrange(permutations):
        perm = 'permutation ' + str(i)
        info(perm)

        afreqs = array(freqs,dtype=float32)
        shuffle(afreqs)

        if missing == -9:
            afreqs = masked_values(afreqs,missing)

        transfreqs = transpose(afreqs)

        for j in range(0,(len(distanceClass))):
            dc = 'distance class ' + str(j)
            info(dc)

            if indice == 'b':
                II, cc = AIDA(transfreqs,DistanceClasses[j])
                MoranII_ci[j].append(II)
                Gearycc_ci[j].append(cc)
            elif indice == 'm':
                II = AIDA(transfreqs,DistanceClasses[j])
                MoranII_ci[j].append(II)
            elif indice == 'g':
                cc = AIDA(transfreqs,DistanceClasses[j])
                Gearycc_ci[j].append(cc)

            #II, cc = AIDA(transfreqs,DistanceClasses[j])
            #MoranII_ci[j].append(II)
            #Gearycc_ci[j].append(cc)


    if indice == 'b' or indice == 'm':
        print """II EMPIRICAL CONFIDENCE INTERVALS, CLASS BY CLASS
(Limits reported are the first and the last value which define the 95, 99 and 99.5% of the empirical area)"""

        i = 0
        for dc in MoranII_ci:
            print "Class #",i

            dc.sort()

            print "95% CI =", (scoreatpercentile(dc, .025)),",",(scoreatpercentile(dc, .975))
            print "99% CI =", (scoreatpercentile(dc, .005)),",",(scoreatpercentile(dc, .995))
            print "99.5% CI =", (scoreatpercentile(dc, .0025)),",",(scoreatpercentile(dc, .9975)),"\n"

            i += 1

    if indice == 'b' or indice == 'g':
        i = 0

        print """cc EMPIRICAL CONFIDENCE INTERVALS, CLASS BY CLASS (Limits reported are the first and the last value which define the 95, 99 and 99.5% of the empirical area)"""


        for dc in Gearycc_ci:

            print "Class #",i

            dc.sort()

            print "95% CI =", (scoreatpercentile(dc, .025)),",",(scoreatpercentile(dc, .975))
            print "99% CI =", (scoreatpercentile(dc, .005)),",",(scoreatpercentile(dc, .995))
            print "99.5% CI =", (scoreatpercentile(dc, .0025)),",",(scoreatpercentile(dc, .9975)),"\n"

            i += 1

elif ci == 's':
    info('begin confidence intervals')

    print "PERMUTATIONS =", p, "\n"

    from random import random, sample
    #from numpy.random import shuffle
    MoranII_ci = []; Gearycc_ci = []

    #from numpy import ones
    Wsum = [sum(Wij) for Wij in DistanceClasses]
    N2 = int(sqrt(min(Wsum)))
    #freqs = sample(freqs,N2)
    Wij = ones((N2, N2))

    for i in xrange(permutations):
        perm = 'permutation ' + str(i)
        info(perm)

        afreqs = sample(freqs,N2)
        afreqs = array(afreqs, dtype=float32)

        if missing == -9:
            afreqs = masked_values(afreqs,missing)

        afreqs2 = (array(afreqs))**2

        if missing == -9:
            afreqs2 = masked_values(afreqs2,missing)

        afreqs2 = afreqs2.sum()

        transfreqs = transpose(afreqs)

        if indice == 'b':
            II, cc = AIDA(transfreqs,Wij,freqs2 = afreqs2,npops = N2)
            MoranII_ci.append(II)
            Gearycc_ci.append(cc)
        elif indice == 'm':
            II = AIDA(transfreqs,Wij,freqs2 = afreqs2,npops = N2)
            MoranII_ci.append(II)
        elif indice == 'g':
            cc = AIDA(transfreqs,Wij,freqs2 = afreqs2,npops = N2)
            Gearycc_ci.append(cc)

    if indice == 'b' or indice == 'm':
        print """II EMPIRICAL CONFIDENCE INTERVALS
(Limits reported are the first and the last value which define the 95, 99 and 99.5% of the empirical area)"""

        MoranII_ci.sort()

        print "95% CI =", (scoreatpercentile(MoranII_ci, .025)),",",(scoreatpercentile(MoranII_ci, .975))
        print "99% CI =", (scoreatpercentile(MoranII_ci, .005)),",",(scoreatpercentile(MoranII_ci, .995))
        print "99.5% CI =", (scoreatpercentile(MoranII_ci, .0025)),",",(scoreatpercentile(MoranII_ci, .9975)),"\n"

    if indice == 'b' or indice == 'g':
        print """cc EMPIRICAL CONFIDENCE INTERVALS
(Limits reported are the first and the last value which define the 95, 99 and 99.5% of the empirical area)"""

        Gearycc_ci.sort()

        print "95% CI =", (scoreatpercentile(Gearycc_ci, .025)),",",(scoreatpercentile(Gearycc_ci, .975))
        print "99% CI =", (scoreatpercentile(Gearycc_ci, .005)),",",(scoreatpercentile(Gearycc_ci, .995))
        print "99.5% CI =", (scoreatpercentile(Gearycc_ci, .0025)),",",(scoreatpercentile(Gearycc_ci, .9975)),"\n"

#gc.collect()
#objgraph.show_most_common_types(limit=20)

sys.stdout = saveout

sys.stdout=sys.__stdout__

#output.close()

info('finish computations')

###------------------------------------Graph time---------------###

"""
from pychart import *
theme.scale_factor = 1
theme.output_format = 'pdf'
theme.output_file = 'MoranII'
theme.reinitialize()

data = zip(distanceClass, MoranII)

xaxis = axis.X(format="/a-60/hL%d", tic_interval = 1, label="Distance Class")
yaxis = axis.Y(label="II")

ar = area.T(x_axis=xaxis, y_axis=yaxis)

plot = line_plot.T(label="II",data=data, tick_mark=tick_mark.star)

ar.add_plot(plot)

ar.draw()

figure(1)
plot(distanceClass,MoranII)
title('II',fontstyle='italic')
xlabel('Distance')
ylabel('II',fontstyle='italic')
savefig(outfile+"_MoranII.png")


xaxis = axis.X(format="/a-60/hL%d", tic_interval = 20, label="Distance Class")
yaxis = axis.Y(tic_interval = 20, label="cc")

figure(2)
plot(distanceClass,Gearycc)
title('cc',fontstyle='italic')
xlabel('Distance Class')
ylabel('cc',fontstyle='italic')
savefig(outfile+"_Gearycc.png")

show()
"""

#for i in range(0,len(MoranII_ci)):
#    print mean(MoranII_ci[i]), mean(Gearycc_ci[i])

#sys.stdout.close( )


