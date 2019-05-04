# -*- coding: utf-8 -*-
import numpy as np
import sys

MIN_THRESHOLD = 3
EXPECTED_FRAGMENT = 29000
QUAL = 40

class Barcode:
    def __init__(self, code):
        self.code = code
        self.Connects = {}
        
    def AddConnect(self, pos, other):
        if other.name in self.Connects:
            self.Connects[other.name].Addpos(pos)
        else:
            self.Connects[other.name] = Connection(pos, other)
            
    def Distribute(self):
        cncts = list(self.Connects.values())
        ct2 = []
        
        kinc = rinc = 0
        
        for cn in cncts:
            kn, rn = cn.AutoFilter()
            kinc += kn
            rinc += rn
            
            if cn.positions:
                ct2.append(cn)
                
        for i in range(len(ct2)):
            for j in range(len(ct2)):
                if i != j:
                    ct2[i].other.CreateConnection(ct2[i], ct2[j].other)
            
        return [kinc, rinc]
        

class Connection:
    def __init__(self, pos, other, allpos = []):
        if allpos:
            self.positions = allpos
        else:
            self.positions = [pos]
        
        self.other = other
        self.weight = 1
        
    def Addpos(self, pos):
        self.positions.append(pos)
        
    def Merge(self, cnct):
        self.weight += 1
        for p in cnct.positions:
            self.positions.append(p)
            
    def AutoFilter(self):
        return self.other.FilterMids(self)
            
class Ark:
    def __init__(self, Connect, other):
        self.positions = Connect.positions.copy()
        
        self.other = other
        self.weight = 1
        
    def Strengthen(self, Connect):
        for p in Connect.positions:
            self.positions.append(p)
        self.weight += 1
            

class Scaffold:
    def __init__(self, snam, length):
        self.name = snam
        self.length = length
        self.Connects = {}
        
    def CreateConnection(self, cnct, other):

        if other.name in self.Connects:
            self.Connects[other.name].Strengthen(cnct)
        else:
            self.Connects[other.name] = Ark(cnct, other)
                
            
    def StripByThreshold(self):
        ndict = {}
        sinc= 0
        ninc = 0
        for oth,cnct in self.Connects.items():
            if cnct.weight >= MIN_THRESHOLD:
                ndict[oth] = cnct
                ninc += cnct.weight
            else:
                 sinc += cnct.weight  
                
        self.Connects = ndict
        
        return [sinc, ninc]
    
    
    def FilterMids(self, cnct):
        ps2 = []
        kinc = rinc = 0
        
        for cn in cnct.positions:
            if cn < EXPECTED_FRAGMENT or (cn > (self.length - EXPECTED_FRAGMENT)):
                ps2.append(cn)
                kinc += 1
            else:
                rinc +=1
        
        cnct.positions = ps2.copy()
        
        return [kinc, rinc]
    
    def GetLRStatus(self, othnm):
        cn = self.Connects[othnm]
        cnm = np.mean(cn.positions)
        ch = 'B'
        
        if self.length > EXPECTED_FRAGMENT:
            if self.length < (EXPECTED_FRAGMENT * 2):
                if cnm > (self.length/2):
                    ch = 'R'
                else:
                    ch = 'L'
            else:
                if (self.length - cnm) >= EXPECTED_FRAGMENT:
                    ch = 'L'
                elif cnm > EXPECTED_FRAGMENT:
                    ch = 'R'
        
        return [ch, cnm]
    
    def CanSplit(self):
        if self.length < EXPECTED_FRAGMENT:
            return False
        else:
            return True

class Genome:
    def __init__(self):
        self.Scaffolds = {}
        self.Barcodes = {}
        self.Candidates = {}
        
    def PrintFasta(self, fname, pref):
        
        cntg = pref + "_candidates.fa"
        
        gn = open(fname, "r")
        go = open(cntg, "w")
        
        swit = 0
        
        for line in gn:
            line = line.rstrip('\r\n')
            
            if line[0] == ">":
                nm = line[1:]
                nm = nm.split(' ')[0]
                if nm in self.Candidates:
                    swit = 1
                else:
                    swit = 0
            
            if swit:
                go.write(line + "\n")
        
        gn.close()
        go.close()
        print("Wrote " + cntg + " successfully!")
     
    def ReadFasta(self, fname):
        glen = 0
        snum = 0
        with open(fname, "r") as gen:
            oldID = ""
            seq = ""
            for line in gen:
                line = line.rstrip('\r\n')
                
                if line[0] == ">":
                    nm = line[1:]
                    nm = nm.split(' ')[0]
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
                
    def ReadPAF(self, fname):
        rinc = linc = 0
        with open(fname, "r") as sfile:
            for line in sfile:
                line = line.rstrip('\r\n')
                
                vls = line.split('\t')
                ql = int(vls[11])
                linc += 1
                
                if ql >= QUAL:
                    nms = vls[0].split('_')
                    sc = vls[5]
                    pos = int(vls[8])
                    coda = nms[-1]
                    
                    if sc in self.Scaffolds:
                        rinc += 1
                        if coda not in self.Barcodes:
                            self.Barcodes[coda] = Barcode(coda)
                        
                        self.Barcodes[coda].AddConnect(pos, self.Scaffolds[sc])
                    
        bl = len(self.Barcodes)
        
        print("Filtered " + str(linc) + " total reads at MAPQ=" + str(QUAL))
        print("Processed " + str(rinc) + " valid alignments")
        print("Found " + str(bl) + " unique barcodes\n")
                        
    def StripScaffolds(self):
        inc1 = inc2 = 0
        for scaff in self.Scaffolds.values():
            inc1a, inc2a = scaff.StripByThreshold()
            inc1 += inc1a
            inc2 += inc2a
            
        print ("This many edges killed: " + str(inc1))
        print ("This many edges kept: " + str(inc2) + "\n")
            
    def DistributeBarcodes(self):
        kinc = rinc = 0
        for bc in self.Barcodes.values():
            kn, rn = bc.Distribute()
            kinc += kn
            rinc += rn
            
        print("Reads kept: " + str(kinc) + 
              "\nReads rejected : " + str(rinc) + "\n")
        
    def AssignCheck(self, ch1, ksi):
        n1 = "nn_" + ksi
        if ch1 != 'B':
            n1 += "_" + ch1
        
        return n1
        
        
    def MirrorDown(self, pref):
        
        tblo = pref + "_table.tsv"
        nno = pref + "_network.tsv"
        ndo = pref + "_nodes.tsv"
        
        ctab = open(tblo, "w")
        nnet = open(nno, "w")
        nnode = open(ndo, "w")
        
        sln = len(self.Scaffolds)
        
        ks = list(self.Scaffolds.keys())
        
        ctab.write("CNTG1\tTYPE1\tMEAN1\tCNTG2\tTYPE2\tMEAN2\tWEIGHT\n")
        nnode.write("NODE\tSIZE\n")
        nnet.write("SOURCE\tWEIGHT\tTARGET\tTYPE\n")
        
        for i in range(sln-1):
            if self.Scaffolds[ks[i]].Connects:
                for j in range(i+1, sln):
                    if ks[j] in self.Scaffolds[ks[i]].Connects:
                        self.Candidates[ks[j]] = 1
                        self.Candidates[ks[i]] = 1
                        
                        wt = self.Scaffolds[ks[i]].Connects[ks[j]].weight
                        ch1, mn1 = self.Scaffolds[ks[i]].GetLRStatus(ks[j])
                        ch2, mn2 = self.Scaffolds[ks[j]].GetLRStatus(ks[i])
                        
                        n1 = self.AssignCheck(ch1, ks[i])
                        n2 = self.AssignCheck(ch2, ks[j])
                            
                        nnet.write(n1 + "\t" + str(wt) + "\t" +
                                   n2 + "\t1\n")

                        ctab.write(ks[i] + "\t" + ch1 + "\t" + str(mn1) +
                                   "\t" + ks[j] + "\t" + ch2 + "\t" + 
                                   str(mn2) + "\t" + str(wt) + "\n")
                        
                
        bl = len(self.Candidates)
        
        for ks in self.Candidates.keys():
            ln = self.Scaffolds[ks].length
            if self.Scaffolds[ks].CanSplit():
                l2 = ln/2
                nnode.write("nn_" + ks + "_L\t" + str(l2) + "\n")
                nnode.write("nn_" + ks + "_R\t" + str(l2) + "\n")
                
                nnet.write("nn_" + ks + "_L\t100\tnn_" + ks + "_R\t2\n")
            else:
                nnode.write("nn_" + ks + "\t" + str(ln) + "\n")
        
        print("Wrote " + tblo + " successfully!")
        print("Wrote " + nno + " successfully!")
        print("Wrote " + ndo + " successfully!\n")
        
        print("Joinable candidate contigs found: " + str(bl) + "\n")
        ctab.close()
        nnet.close()
        nnode.close()
                        

def HelpMsg():
    print ("\nNanochrome - building barcode linkage candidates\n")
    print ("ARG1: <genome.fa> A Genome file")
    print ("ARG2: <alns.paf> Mapped 10X library")
    print ("ARG3: <integer> Fragment Length (x10 library)")
    print ("ARG4: <prefix> Unique prefix for outputs\n")
    print ("Output:\n")
    print ("<prefix>_candidates.fa: A set of scaffold candidates")
    print ("<prefix>_table.tsv: Potential connection table")
    print ("<prefix>_network.tsv: A network file for visualisation")
    print ("<prefix>_nodes.tsv: A node description file for visualisation\n")

def main(argv):
    if len(argv) >= 4:
        try:
            global EXPECTED_FRAGMENT
            EXPECTED_FRAGMENT = int(argv[2])
            print("Reading Genome...\n")
            the_genome = Genome()
            the_genome.ReadFasta(argv[0])
            print("Reading Alignments...\n")
            the_genome.ReadPAF(argv[1])
            print("Filtering Reads & Edges...\n")
            the_genome.DistributeBarcodes()
            the_genome.StripScaffolds()
            print("Generating Barcode Graph..\n")
            the_genome.MirrorDown(argv[3])
            print("Writing Candidates...\n")
            the_genome.PrintFasta(argv[0], argv[3])
            #sys.exit(2)
        except:
            print("Error: ",sys.exc_info()[0]," <- this happened.")
            HelpMsg()
        
    else:
        HelpMsg()
        print("You need to specify the four positional arguments\n")
        #sys.exit(2)
 
 
if __name__ == '__main__':
    main(sys.argv[1:])