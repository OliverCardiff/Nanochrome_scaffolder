# -*- coding: utf-8 -*-

import sys
import numpy as np

DEGREE_FILTER = 5
QUAL = 25

class Join:
    def __init__(self, other, weight):
        self.other = other
        self.weight = weight
        self.nano_confs = 0
        self.orientation = 'F'
        self.gap = 500
    
class Scaffold:
    def __init__(self, name, length):
        self.name = name
        self.length = length
        self.bestpath = 0
        self.seq = ""
        
        self.Left = {}
        self.Right = {}
        self.Both = {}
        
    def AddJoin(self, other, weight, mode):
        jn = Join(other, weight)
        
        if mode == 'L':
            self.Left[other.name] = jn
        elif mode == 'R':
            self.Right[other.name] = jn
        else:
            self.Both[other.name] = jn
            
    def CanRecieveNano(self, aln, rlength, other):
        
        mds = {'L':0,'R':0,'B':0}
        resids = {'L':0, 'R':0, 'B':0}
        
        if other.name in self.Left:
            mds['L'] += 1
        if other.name in self.Right:
            mds['R'] += 1
        if other.name in self.Both:
            mds['B'] += 1
            
        if mds['L'] or mds['B']:
            mlen = rlen = 0
            if aln.rev:
                mlen = aln.refend
                rlen = rlength - aln.qend
            else:
                mlen = aln.refstart
                rlen = aln.qstart
            
            resids['L'] = rlen - mlen
            if rlen > mlen:
                mds['L'] += 1
                
        if mds['R'] or mds['B']:
            mlen = rlen = 0
            if aln.rev:
                mlen = self.length - aln.refstart
                rlen = rlength - aln.qend
            else:
                mlen = self.length - aln.refend
                rlen = aln.qend
            
            resids['R'] = rlen - mlen
            if rlen > mlen:
                mds['R'] += 1
                
        if mds['B']:
            resid = max([resids['R'], resids['L']])
            if resid > 0:
                mds['B'] += 1
                
        if max([mds.values()]) == 2:
            return True
        
        return False
    
    def __stuncf(self, dct):
        ks = list(dct.keys())
        for k in ks:
            if not dct[k].nano_confs:
                del dct[k]
                
    def StripUnconf(self):
        self.__stuncf(self.Left)
        self.__stuncf(self.Right)
        self.__stuncf(self.Both)
            
    def DegreeFilter(self):
        hval = 0
        def FilterOut(dct):
            for jns in dct.values():
                jns.other.KillOther(self.name)
                
        if len(self.Left) > DEGREE_FILTER:
            hval += len(self.Left)
            FilterOut(self.Left)
            self.Left = {}
            
        if len(self.Right) > DEGREE_FILTER:
            hval += len(self.Right)
            FilterOut(self.Right)
            self.Right = {}
            
        if len(self.Both) > DEGREE_FILTER:
            hval += len(self.Both)
            FilterOut(self.Both)
            self.Both = {}
            
        return hval
            
            
    def KillOther(self, name):
        if name in self.Right:
            del self.Right[name]
        
        if name in self.Left:
            del self.Left[name]
            
        if name in self.Both:
            del self.Both[name]

class Alignment:
    def __init__(self, qst, qed, rst, red, rev):
        self.qstart = qst
        self.qend = qed
        self.refstart = rst
        self.refend = red
        self.rev = rev
        
class Read:
    def __init__(self, coda, ln):
        self.coda = coda
        self.alignments = {}
        self.length = ln
    
    def AddAlign(self, qst, qed, rst, red, rev, othname):
        aln = Alignment(qst, qed, rst, red, rev)
        self.alignments[othname] = aln

class Graph:
    def __init__(self):
        self.Scaffolds = {}
        self.Reads = {}
        
    def StripUnconfirmed(self):
        for r in self.Reads.values():
            r.StripUnconf()
        
    def PreFilterEdges(self):
        hkill = 0
        ekill = 0
        for scs in self.Scaffolds.values():
            res = scs.DegreeFilter()
            hkill += 1 if res > 0 else 0
            ekill += res
        
        print("Killed: " + str(hkill) + " hubs over degree of " +
              str(DEGREE_FILTER))
        print("Removing a total of " + ekill + " edges\n")
        
    def ScanGenome(self, fname):
        glen = 0
        snum = 0
        with open(fname, "r") as gen:
            oldID = ""
            seq = ""
            for line in gen:
                line = line.rstrip('\r\n')
                
                if line[0] == ">":
                    nm = line[1:]
                    if len(oldID) > 0 and len(seq) > 150:
                        self.Scaffolds[oldID] = Scaffold(oldID, len(seq))
                        glen += len(seq)
                        snum += 1
                        seq = ""
                        
                    oldID = nm
                else:
                    seq += line
            else:
                self.Scaffolds[oldID] = Scaffold(oldID, len(seq))
                glen += len(seq)
                snum += 1
                
        print("Read Genome with " + str(glen) + " bases")
        print("Found " + str(snum) + " separate contigs\n")
    
    def ReadTable(self, tnam):
        
        def GetInc(t1, t2, tst):
            rv = 0
            if t1 == tst:
                rv += 1
                
            if t2 == tst:
                rv += 1
                
            return rv
        
        jinc = linc = rinc = binc = 0
        
        with open(tnam, "r") as tblin:
            next(tblin)
            for line in tblin:
                line.rstrip('\r\n')
                lns = line.split('\t')
                
                wt = lns[4]
                t1 = lns[1]
                t2 = lns[3]
                
                jinc += 2
                rinc += GetInc(t1,t2, 'R')
                linc += GetInc(t1,t2, 'L')
                binc += GetInc(t1,t2, 'B')
                
                if lns[0] in self.Scaffolds and lns[2] in self.Scaffolds:
                    sc1 = self.Scaffolds[lns[0]]
                    sc2 = self.Scaffolds[lns[2]]
                    sc1.AddJoin(sc2, wt, t1)
                    sc2.AddJoin(sc1, wt, t2)

        print ("Loaded " + jinc + " total joins")
        print ("With " + linc + " left arks")
        print ("With " + rinc + " right arks")
        print ("With " + binc + " ambiguous arks\n")
        
    def ReadPAF(self, fname):
        rinc = linc = revs = uniqs = 0
        currRead = None
        with open(fname, "r") as sfile:
            for line in sfile:
                line = line.rstrip('\r\n')
                
                if line[0] != '@':
                    vls = line.split('\t')
                    
                    coda = vls[0]
                    clen = int(vls[1])
                    ast = int(vls[2])
                    aed = int(vls[3])
                    rev = vls[4]
                    sc = vls[5]
                    tst = int(vls[7])
                    ted = int(vls[8])                    
                    ql = int(vls[11])

                    linc += 1
                    
                    bit = 1 if rev == '+' else 0
                    
                    if sc in self.Scaffolds and 255 != ql >= QUAL:
                        rinc += 1
                        revs += 1 if bit else 0
                        if not currRead:
                            currRead = Read(coda, clen)
                        if coda != currRead.coda:
                            currRead = Read(coda, clen)
                            self.Reads[coda] = currRead
                            
                        currRead.AddAlign(ast, aed, tst, ted, bit, sc)
                            
                        

        uniqs = len(self.Reads)        
        print("Filtered " + str(linc) + " total aligns at MAPQ=" + str(QUAL))
        print("Processed " + str(rinc) + " valid alignments")
        print("..of which " + str(revs) +  " were reversed")
        print("From a total of " + str(uniqs) + " unique reads\n")
        
    def DistributeReads(self):
        for rd in self.Reads.values():
            a_len = len(rd.alignments)
            rln = rd.length
            if a_len > 1:
                ks = list(rd.alignments.keys())
                als = list(rd.alignments.values())
                for i in range(a_len-1):
                    for j in range(i+1,a_len):
                        sc1 = self.Scaffolds[ks[i]]
                        sc2 = self.Scaffolds[ks[j]]
                        if sc1.CanRecieveNano(als[i], rln, sc2):
                            if sc2.CanRecieveNano(als[j], rln, sc1):
                                self.AttemptToReconcile(sc1,sc2,
                                                        als[i],als[j],rln)
                                
                            
                            
    def AttemptToReconcile(self, scaff1, scaff2, aln1, aln2, rlen):
        
        return True
    
    def StoreGenome(self, fname):
        print("Reloading candidate sequences...\n")
        with open(fname, "r") as gen:
            oldID = ""
            seq = ""
            for line in gen:
                line = line.rstrip('\r\n')
                
                if line[0] == ">":
                    nm = line[1:]
                    if len(oldID) > 0 and oldID in self.Scaffolds:
                        self.Scaffolds[oldID].seq = seq
                        seq = ""
                        
                    oldID = nm
                else:
                    seq += line
            else:
                if oldID in self.Scaffolds:
                    self.Scaffolds[oldID].seq = seq
                
        
def HelpMsg():
    print ("NanoChrome - Confirming 10x with nanopore\n\n")
    print ("ARG1: <genome.fa> A Genome file\n")
    print ("ARG2: <nc_table.tsv> output from chrome_candidates.py")
    print ("ARG3: <alns.sam> Mapped Longread library\n")
    print ("ARG4: <prefix> Unique prefix for outputs\n\n")
    print ("Output:\n")
    print ("<prefix>_nc_scaffolded.fa: Final Scaffolded Genome\n")
    print ("<prefix>_network_final.tsv: A network file for visualisation\n")
    print ("<prefix>_nodes_final.tsv: A node description file for visualisation\n\n")

def main(argv):

    # make sure there are at least three arguments
    if len(argv) >= 4:
        try:
            the_graph = Graph()
            print("Scanning the genome...\n")
            the_graph.ScanGenome(argv[0])
            print("Reading the edge graph...\n")
            the_graph.ReadTable(argv[1])
            print("Pre-filtering edge graph...\n")
            the_graph.PreFilterEdges()
            print("Reading in longread PAF...\n")
            the_graph.ReadPAF(argv[2])
            print("Distributing Reads to graph...\n")
            the_graph.DistributeReads()
            print("Agent navigating graph...\n")
            the_graph.StripUnconfirmed()
            print("Applying Optimal Scaff Paths\n")
            the_graph.StoreGenome(argv[0])
            print("Writing Final Assembly\n")
            
            
        except:
            print("Error: ",sys.exc_info()[0]," <- this happened.")
        finally:
           HelpMsg() 
        
    else:
        HelpMsg()
        print("You need to specify the four positional arguments\n")
        sys.exit(2)
 
 
if __name__ == '__main__':
    main(sys.argv[1:])