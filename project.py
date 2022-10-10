# HEDGES Error-Correcting Code for DNA Storage Corrects Indels and Allows Sequence Constraints
# William H. Press, John A. Hawkins, Stephen Knox Jones Jr, Jeffrey M. Schaub, Ilya J. Finkelstein
# submitted to Proceedings of the National Academy of Sciences
#
# Demonstration driver and installation validation program
#
# This file demonstrates the use of the separately compiled modules NRpyDNAcode and NRpyRS.
# (Those modules are compiled in C++ using Numerical Recipes.  See separate documentation.)
#
# We encode a specified number of packets from known plaintext, create DNA errors, then
# decode the DNA and compare.
import numpy.version
from numpy import *
import NRpyDNAcode as code
import NRpyRS as RS
import sys
import pathlib
import os
import datetime

#mydir = "D:\\Dropbox\\Projects\\DNAcode\\CodeAsPublished"
#os.chdir(mydir)
#print help(code) # uncomment to see NRpyDNAcode help output
coderates = array([NaN, 0.75, 0.6, 0.5, 1./3., 0.25, 1./6.]) # table of coderates 1..6

# user-settable parameters for this test
coderatecode = 3 # test this coderate in coderates table above
npackets = 20 # number of packets (of 255 strands each) to generate and test
totstrandlen = 300 # total length of DNA strand
strandIDbytes = 2 # ID bytes each strand for packet and sequence number
strandrunoutbytes = 2 # confirming bytes end of each strand (see paper)
hlimit = 1000000 # maximum size of decode heap, see paper
leftprimer = "TCGAAGTCAGCGTGTATTGTATG"
rightprimer = "TAGTGAGTGCGATTAAGCGTGTT" # for direct right appending (no revcomp)

# this test generates substitution, deletion, and insertion errors
# sub,del,ins rates to simulate (as multiple of our observed values):

totalStrandLenOption = int(input("Press 0 for default [300] total strand length of the DNA (including left and right primers)"
                                 ", 1 for custom length: "))
if totalStrandLenOption == 1:
    strandLenCandidate = int(input("Total strand length of the DNA (must be more than 46): "))
    if strandLenCandidate > 46:
        totstrandlen = strandLenCandidate

outputPathOption = int(input("Press 0 for default [stdout] output path, 1 for custom path: "))
if outputPathOption == 1:
    outputPathCandidate = input("output path")

codeRateOption = int(input("Press 0 for default [0.5] code rate, 1 for custom code rate: "))
if codeRateOption == 1:
    coderatecode = int(input("Press 1 for 0.75, 2 for 0.6, 3 for 0.5, 4 for 0.33, 5 for 0.25, 6 for 0.166: "))


ratesOption = int(input("Press 0 for default [s:0.0238, d:0.0082, i:0.0039] Substitution/Deletion/Insertion rates, 1 for custom rates: "))
if ratesOption == 0:
    (srate, drate, irate) = 1.5 * array([0.0238, 0.0082, 0.0039])
else:
    srate = float(input("Substitution Rate: "))
    drate = float(input("Deletion Rate: "))
    irate = float(input("Insertion Rate: "))

# set parameters for DNA constraints (normally not changed, except for no constraint)
max_hpoly_run = 4 # max homopolymer length allowed (0 for no constraint)
GC_window = 12 # window for GC count (0 for no constraint)
max_GC = 8 # max GC allowed in window (0 for no constraint)
min_GC = GC_window - max_GC

# not normally user settable because assumed by Reed-Solomon outer code:
strandsperpacket = 255  # for RS(255,32)
strandsperpacketcheck = 32  # for RS(255,32)

# compute some derived parameters and set parameters in NRpyDNAcode module
leftlen = len(leftprimer)
rightlen = len(rightprimer)
strandlen = totstrandlen - leftlen - rightlen
strandsperpacketmessage = strandsperpacket - strandsperpacketcheck
(NSALT, MAXSEQ, NSTAK, HLIMIT) = code.getparams() # get settable code parameters
code.setparams(8*strandIDbytes, MAXSEQ, NSTAK, hlimit) # change NSALT and HLIMIT
bytesperstrand = int(strandlen*coderates[coderatecode]/4.)    
messbytesperstrand = bytesperstrand - strandIDbytes - strandrunoutbytes # payload bytes per strand
messbytesperpacket = strandsperpacket * messbytesperstrand # payload bytes per packet of 255 strands
code.setcoderate(coderatecode, leftprimer, rightprimer) # set code rate with left and right primers
code.setdnaconstraints(GC_window, max_GC, min_GC, max_hpoly_run) # set DNA constraints (see paper)

# define a source of plaintext bytes, either random or The Wizard of Oz in Esperanto
dataOption = int(input("Press 0 for default [WizardOfOzInEsperanto.txt] input data file, 1 for custom data file, 2 for binary file: "))
if dataOption == 1 or dataOption == 2:
    customDataInputFile = input('Data file name: ')
else:
    customDataInputFile = "WizardOfOzInEsperanto.txt"
print('Using file: {} as data input'.format(customDataInputFile))

fileoffset = 0
if dataOption == 2:
    with open(customDataInputFile, 'rb') as myfile: filetext = myfile.read()
    filebytes = array([c for c in filetext])
else:
    with open(customDataInputFile, 'r') as myfile: filetext = myfile.read()
    filebytes = array([ord(c) for c in filetext]).astype(uint8)
wizlen = len(filebytes)
def getdatafile(n): # return next n chars from wiztext
    global fileoffset, filelen
    if fileoffset + n > wizlen : fileoffset = 0
    bytes = filebytes[fileoffset:fileoffset+n]
    fileoffset += n
    return bytes

# functions to create sequential packets from the plaintext source, and R-S protect them
def createmesspacket(packno) : # packno in range 0..255 with value 2 for strandIDbytes
    packet = zeros([strandsperpacket,bytesperstrand], dtype=uint8)
    plaintext = zeros(strandsperpacketmessage*messbytesperstrand, dtype=uint8)
    for i in range(strandsperpacket):
        packet[i, 0] = packno # note assumes value 2 for strandIDbytes
        packet[i, 1] = i
        if i < strandsperpacketmessage :
            ptext = getdatafile(messbytesperstrand)
            packet[i, strandIDbytes:strandIDbytes+messbytesperstrand] = ptext
            plaintext[i*messbytesperstrand:(i+1)*messbytesperstrand] = ptext
    return packet, plaintext

def protectmesspacket(packetin):  # fills in the RS check strands
    packet = packetin.copy()
    regin = zeros(strandsperpacket, dtype=uint8)
    for j in range(messbytesperstrand):
        for i in range(strandsperpacket):
            regin[i] = packet[i, ((j+i) % messbytesperstrand)+strandIDbytes]
        regout = RS.rsencode(regin)
        for i in range(strandsperpacket):
            packet[i, ((j+i) % messbytesperstrand)+strandIDbytes] = regout[i]
    return packet

# functions to encode a packet to DNA strands, and decode DNA strands to a packet
def messtodna(mpacket):
    # HEDGES encode a message packet into strands of DNA
    filler = array([0, 2, 1, 3, 0, 3, 2, 1, 2, 0, 3, 1, 3, 1, 2, 0, 2, 3, 1, 0, 3, 2, 1, 0, 1, 3], dtype=uint8)
    dpacket = zeros([strandsperpacket, totstrandlen], dtype=uint8)
    for i in range(strandsperpacket):
        dna = code.encode(mpacket[i, :])
        if len(dna) < totstrandlen: # need filler after message and before right primer
            dnaleft = dna[:-rightlen]
            dnaright = dna[-rightlen:]
            dna = concatenate((dnaleft, filler[:totstrandlen-len(dna)],dnaright))
            #n.b. this can violate the output constraints (very slightly at end of strand)
        dpacket[i, :len(dna)] = dna
    return dpacket

def dnatomess(dnapacket):
    # HEDGES decode strands of DNA (assumed ordered by packet and ID number) to a packet
    baddecodes = 0
    erasures = 0
    mpacket = zeros([strandsperpacket, bytesperstrand], dtype=uint8)
    epacket = ones([strandsperpacket, bytesperstrand], dtype=uint8)  # everything starts as an erasure
    for i in range(strandsperpacket):
        (errcode, mess, _, _, _, _) = code.decode(dnapacket[i, :], 8*bytesperstrand)
        if errcode > 0:
            baddecodes += 1
            erasures += max(0,messbytesperstrand-len(mess))
        lenmin = min(len(mess),bytesperstrand)
        mpacket[i, :lenmin] = mess[:lenmin]
        epacket[i, :lenmin] = 0
    return mpacket, epacket, baddecodes, erasures

# functions to R-S correct a packet and extract its payload to an array of bytes
def correctmesspacket(packetin, epacket):
    # error correction of the outer RS code from a HEDGES decoded packet and erasure mask
    packet = packetin.copy()
    regin = zeros(strandsperpacket, dtype=uint8)
    erase = zeros(strandsperpacket, dtype=uint8)
    tot_detect = 0
    tot_uncorrect = 0
    max_detect = 0
    max_uncorrect = 0
    toterrcodes = 0
    for j in range(messbytesperstrand):
        for i in range(strandsperpacket):
            regin[i] = packet[i, ((j+i) % messbytesperstrand)+strandIDbytes]
            erase[i] = epacket[i, ((j+i) % messbytesperstrand)+strandIDbytes]
        locations = array(argwhere(erase), dtype=int32)
        (decoded, errs_detected, errs_corrected, errcode, ok) = RS.rsdecode(regin, locations)
        tot_detect += errs_detected
        tot_uncorrect += max(0, (errs_detected-errs_corrected))
        max_detect = max(max_detect, errs_detected)
        max_uncorrect = max(max_uncorrect, max(0, (errs_detected-errs_corrected)))
        toterrcodes += (0 if errcode == 0 else 1)
        for i in range(strandsperpacket):
            packet[i, ((j+i) % messbytesperstrand)+strandIDbytes] = decoded[i]
    return packet, tot_detect, tot_uncorrect, max_detect, max_uncorrect, toterrcodes

def extractplaintext(cpacket):
    # extract plaintext from a corrected packet
    plaintext = zeros(strandsperpacketmessage*messbytesperstrand, dtype=uint8)
    for i in range(strandsperpacketmessage):
        plaintext[i*messbytesperstrand:(i+1)*messbytesperstrand] = (
            cpacket[i, strandIDbytes:strandIDbytes+messbytesperstrand])
    return plaintext

# function to create errors in a bag (or packet) of DNA strands
def createerrors(dnabag,srate,drate,irate) :
    # for testing: create errors in a bag of strands
    (nrows,ncols) = dnabag.shape
    newbag = zeros([nrows, ncols],dtype=uint8)
    for i in range(nrows):
        dna = code.createerrors(dnabag[i, :], srate, drate, irate)
        lenmin = min(len(dna), ncols)
        newbag[i, :lenmin] = dna[:lenmin]
    return newbag

if outputPathOption == 1:
    sys.stdout = open(outputPathCandidate, "w")
    print('Using file: {} as output path'.format(outputPathCandidate))
else:
    print('Using default stdout as output path')
## DO THE TEST
print("DNA mapping used in the project [0: A, 1: C, 2: G, 3: T]")
print("------------------------------------------")
print("Log file for the current run can be found under logs directory with a timestamp of the current time")
print("the log file includes for each packet: packet number, message packet, plain message, reed solomon packet, dna packet.")
print("------------------------------------------")
print("for each packet, these statistics are shown in two groups:")
print("Packet number: (1.1, 1.2, 1.3, 1.4) (2.1, 2.2, 2.3, 2.4) packet OK/NOT")
print("------------------------------------------")
print("1.1 HEDGES decode failures")
print("1.2 HEDGES bytes thus declared as erasures")
print("1.3 R-S total errors detected in packet")
print("1.4 max errors detected in a single decode")
print("------------------------------------------")
print("2.1 R-S reported as initially-uncorrected-but-recoverable total")
print("2.2 same, but max in single decode")
print("2.3 R-S total error codes; if zero, then R-S corrected all errors")
print("2.4 Actual number of byte errors when compared to known plaintext input")
print("------------------------------------------")
badpackets = 0
Totalbads = zeros(8, dtype=int)

# Create a logging directory
LogsPath = str(pathlib.Path().resolve()) + "/logs"
pathlib.Path(LogsPath).mkdir(parents=True, exist_ok=True)

# Create a log for the current run
currentTime = str(datetime.datetime.now())
currentTime = currentTime.replace(':', '_')
file = open(LogsPath + "/" + currentTime + ".txt",'a')

file.write('HEDGES error-correcting code for DNA storage corrects indels and allows sequence constraints\n')
numpy.set_printoptions(threshold=sys.maxsize)
for ipacket in range(npackets) :
    file.write('--------------------------------------------------------------------------------------------\n')
    file.write("Packet number: %d \n" % ipacket)
    file.write("--------------------\n")
    # encode
    messpack, messplain = createmesspacket(ipacket)  # plaintext to message packet
    file.write("messpack: \n")
    file.write("--------------------\n")
    file.writelines(str(messpack) + "\n")
    file.write("--------------------\n")
    file.write("messplain: \n")
    file.write("--------------------\n")
    file.writelines(str(messplain)+ "\n")
    file.write("--------------------\n")
    rspack = protectmesspacket(messpack)  # Reed-Solomon protect the packet
    file.write("rspack: \n")
    file.write("--------------------\n")
    file.writelines(str(rspack)+ "\n")
    file.write("--------------------\n")
    dnapack = messtodna(rspack)  # encode to strands of DNA containing payload messplain
    file.write("dnapack: \n")
    file.write("--------------------\n")
    file.writelines(str(dnapack)+ "\n")
    file.write("--------------------\n")
    # simulate errors in DNA synthesis and sequencing
    obspack = createerrors(dnapack, srate, drate, irate)

    # decode
    (dpacket, epacket, baddecodes, erasures) = dnatomess(obspack) # decode the strands
    (cpacket, tot_detect, tot_uncorrect, max_detect, max_uncorrect, toterrcodes) = correctmesspacket(dpacket, epacket)
    
    # check against ground truth
    messcheck = extractplaintext(cpacket)
    badbytes = count_nonzero(messplain-messcheck)
    
    # print summary line
    Totalbads += array([baddecodes, erasures, tot_detect, max_detect, tot_uncorrect, max_uncorrect, toterrcodes, badbytes])
    print(("%3d: (%3d %3d %3d %3d) (%3d %3d %3d %3d)" % (ipacket, baddecodes, erasures,
        tot_detect, max_detect, tot_uncorrect, max_uncorrect, toterrcodes, badbytes)), end=' ')
    print("packet OK" if badbytes == 0 else "packet NOT ok")
    if badbytes: badpackets += 1
print("all packets OK" if not badpackets else "some packets had errors!")
print("TOT: (%4d %4d %4d %4d) (%4d %4d %4d %4d)" % tuple(Totalbads))
