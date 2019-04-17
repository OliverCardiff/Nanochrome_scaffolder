# -*- coding: utf-8 -*-

import sys
import pandas as pd

QUANTILE_THRESH = 0.97
CONF_WEIGHT = 30
DEGREE_FILTER = 5
QUAL = 25
AUTO_GAP = 2000
TRIM_LIMIT = 500

class Join:
    def __init__(self, other, weight):
        self.other = other
        self.weight = weight
        self.nano_confs = 0
        self.auto_confs = 0
        self.orientation = 'F'
        self.gap = AUTO_GAP
        self.trim_amnt = 0
        self.best = False
    
    def __lt__(self,other):
        return self.weight < other.weight
    
class Scaffold:
    def __init__(self, name, length):
        self.name = name
        self.length = length
        self.bestpath = 0
        self.seq = ""
        
        self.Left = {}
        self.Right = {}
        self.Both = {}
    
    def SwitchDict(self, m1, m2, othnm):
        def SendJoin(jn, oth, m):
            if m == 'L':
                self.Left[oth] = jn
            elif m == 'R':
                self.Right[oth] = jn
            else:
                self.Both[oth] = jn
                
        if m1 != m2:
            if m1 == 'B':
                SendJoin(self.Both[othnm], othnm, m2)
            elif m1 == 'L':
                SendJoin(self.Left[othnm], othnm, m2)
            elif m1 == 'R':
                SendJoin(self.Right[othnm], othnm, m2)
    
    def AnchorNano(self, othnm, old_mode, new_mode, trim, gap, rev):
        def UpdateSwitch(jn, trim, gap, rev):
            jn.trim_amnt = trim
            jn.gap = gap
            jn.nano_confs = 1
            if rev:
                jn.orientation = 'R'
                
        self.SwitchDict(old_mode, new_mode, othnm)
        
        if new_mode == 'L':
            UpdateSwitch(self.Left[othnm], trim, gap, rev)
            
        elif new_mode == 'R':
            UpdateSwitch(self.Right[othnm], trim, gap, rev)
        
    
    def HasFixedMatch(self, scfnm):
        if scfnm in self.Left:
            return 'L'
        elif scfnm in self.Right:
            return 'R'
        else:
            return 'N'
             
    def RunAutoConf(self):
        def CheckEntry(vl, ch1, ch2):
            if vl.weight > CONF_WEIGHT:
                res = vl.other.HasFixedMatch(self.name)
                if res == 'L':
                    vl.orientation = ch1
                    vl.auto_confs = 1
                    return 1
                if res == 'R':
                    vl.orientation = ch2
                    vl.auto_confs = 1
                    return 1
            return 0
        
        cnf = 0
        allc = 0
        for ks,vl in self.Left.items():
            allc += 1
            cnf += CheckEntry(vl, 'R', 'F')
            
        for ks,vl in self.Right.items():
            allc += 1
            cnf += CheckEntry(vl, 'F', 'R')
        
        allc += len(self.Both)
        
        return [allc, cnf]
        
    def AddJoin(self, other, weight, mode):
        jn = Join(other, weight)
        
        if mode == 'L':
            self.Left[other.name] = jn
            self.Left[other.name].mode = 'L'
        elif mode == 'R':
            self.Right[other.name] = jn
            self.Right[other.name].mode = 'R'
        else:
            self.Both[other.name] = jn
            self.Both[other.name].mode = 'B'
            
    def GetTail(self, aln, mode):
        res = 0
        
        if mode == 'L':
            if aln.rev:
                res = aln.refend
            else:
                res = aln.refstart
        elif mode == 'R':
            if aln.rev:
                res = self.length - aln.refstart
            else:
                res = self.length - aln.refend
        
        return res
            
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
                
        if max(list(mds.values())) > 1:
            return True
        
        return False
    
    def GetJoinStatus(self, othnm):
        if othnm in self.Left:
            return 'L'
        elif othnm in self.Right:
            return 'R'
        elif othnm in self.Both:
            return 'B'
        else:
            return 'N'
    
    def CountEdges(self):
        return len(self.Left) + len(self.Right) + len(self.Both)
    
    def __stuncf(self, dct):
        removed = 0
        ks = list(dct.keys())
        for k in ks:
            if not dct[k].nano_confs and not dct[k].auto_confs:
                del dct[k]
                removed += 1
        return removed
                
    def StripUnconf(self):
        rem = 0
        nrem = self.CountEdges()
        rem += self.__stuncf(self.Left)
        rem += self.__stuncf(self.Right)
        rem += self.__stuncf(self.Both)
        nrem = nrem - rem
        return [rem, nrem]
            
    def DegreeFilter(self):
        hval = 0
        def FilterOut(dct):
            for jns in dct.values():
                jns.other.KillOther(self.name)
        
        def SplitDict(dct):
            keep = {}
            chuck = {}
            
            stupes = sorted(dct.items(), key=lambda kv: kv[1], reverse=True)
            
            for i in range(len(stupes)):
                if i < DEGREE_FILTER:
                    keep[stupes[i][0]] = stupes[i][1]
                else:
                    chuck[stupes[i][0]] = stupes[i][1]
            
            return [keep, chuck]
        
        def FilterDict(dct):
            replacement = dct
            hval = 0
            ln = len(dct)
            
            if ln > DEGREE_FILTER:
                if ln > (DEGREE_FILTER * 10):
                    hval = ln
                    FilterOut(dct)
                    replacement = {}
                else:
                    hval = ln - DEGREE_FILTER
                    k, c = SplitDict(dct)
                    replacement = k
                    FilterOut(c)
                    
            return replacement, hval
        
        ret, hv = FilterDict(self.Left)
        
        self.Left = ret
        hval += hv
        
        ret, hv = FilterDict(self.Right)
        
        self.Right = ret
        hval += hv
        
        ret, hv = FilterDict(self.Both)
        
        self.Both = ret
        hval += hv
            
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
        self.matched = 0
    
    def AddAlign(self, qst, qed, rst, red, rev, othname):
        aln = Alignment(qst, qed, rst, red, rev)
        self.alignments[othname] = aln

class Graph:
    def __init__(self):
        self.Scaffolds = {}
        self.Reads = {}
        
    def StripUnconfirmed(self):
        rem  = nrem = scs = 0
        for sc in self.Scaffolds.values():
            rem1, nrem1 = sc.StripUnconf()
            rem += rem1
            nrem += nrem1
            scs += 1 if nrem1 else 0
            
        hlf1 = int(rem/2)
        hlf2 = int(nrem/2)
        print("Stripped " + str(hlf1) + " unconfirmed edges from graph")
        print("..leaving " + str(hlf2) + " edges in graph")
        print("..shared by " + str(scs) + " scaffolds\n")
        
    def PreConfirmEdges(self):
        allc = cnf = 0
        
        for scs in self.Scaffolds.values():
            allc1, cnf1 = scs.RunAutoConf()
            allc += allc1
            cnf += cnf1
        
        hfl = int(cnf/2)
        print("Scanned " + str(allc) + " outgoing edges")
        print("Preconfirmed " + str(hfl) + " long frag joins with " +
              str(CONF_WEIGHT) + " barcode supports\n")
        
    def PreFilterEdges(self):
        hkill = 0
        ekill = 0
        for scs in self.Scaffolds.values():
            res = scs.DegreeFilter()
            hkill += 1 if res > 0 else 0
            ekill += res
        
        print("Shrunk: " + str(hkill) + " hubs over degree of " +
              str(DEGREE_FILTER))
        print("Removing a total of " + str(ekill) + " edges\n")
        
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
        
        tbl = pd.read_csv(tnam, sep="\t")
        print("Actual join count: " + str(tbl.shape[0]))
        CONF_WEIGHT = tbl.WEIGHT.quantile(QUANTILE_THRESH)
        print(str(QUANTILE_THRESH) + " quantile filter; weight threshold: " +
              str(CONF_WEIGHT))
        tbl = tbl.loc[tbl.WEIGHT > CONF_WEIGHT]
        print("Post filter count to load: " + str(tbl.shape[0]) + "\n")
        
        for i in range(len(tbl)):
            row = tbl.iloc[i]
            
            wt = row.WEIGHT
            t1 = row.TYPE1
            t2 = row.TYPE2
            sc1nm = str(row.CNTG1)
            sc2nm = str(row.CNTG2)
            
            jinc += 2
            rinc += GetInc(t1,t2, 'R')
            linc += GetInc(t1,t2, 'L')
            binc += GetInc(t1,t2, 'B')
            
            if sc1nm in self.Scaffolds and sc2nm in self.Scaffolds:
                sc1 = self.Scaffolds[sc1nm]
                sc2 = self.Scaffolds[sc2nm]
                sc1.AddJoin(sc2, wt, t1)
                sc2.AddJoin(sc1, wt, t2)

        print ("Loaded " + str(jinc) + " partial joins")
        print ("With " + str(linc) + " left arks")
        print ("With " + str(rinc) + " right arks")
        print ("With " + str(binc) + " ambiguous arks\n")
        
    def ReadPAF(self, fname):
        rinc = linc = revs = uniqs = confirmed = removed = 0
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
                    
                    if rev == '-':
                        tmp = tst
                        tst = ted
                        ted = tmp

                    linc += 1
                    
                    bit = 1 if rev == '-' else 0
                    
                    if linc % 10000 == 0:
                        confirmed += self.DistributeReads()
                        removed += self.FilterReads()
                    
                    if sc in self.Scaffolds and 255 != ql >= QUAL:
                        rinc += 1
                        revs += 1 if bit else 0
                        if not currRead:
                            currRead = Read(coda, clen)
                        if coda != currRead.coda:
                            currRead = Read(coda, clen)
                            self.Reads[coda] = currRead
                            
                        currRead.AddAlign(ast, aed, tst, ted, bit, sc)
                            
                        
        confirmed += self.DistributeReads()
        removed += self.FilterReads()
        uniqs = len(self.Reads)      
        totalr = uniqs + removed
        
        print("Filtered " + str(linc) + " total aligns at MAPQ=" + str(QUAL))
        print("Processed " + str(rinc) + " valid alignments")
        print("..of which " + str(revs) +  " were reversed")
        print("From a total of " + str(totalr) + " unique reads")
        print("..this many: " + str(removed) + ", were unsuitable")
        print("..this many: " + str(uniqs) + ", were useful!\n")
        print("Confirmed " + str(confirmed) + " chromium edges with longreads\n")
        
    def FilterReads(self):
        rem = 0
        ks = list(self.Reads.keys())
        for k in ks:
            if not self.Reads[k].matched or len(self.Reads[k].alignments) < 2:
                del self.Reads[k]
                rem += 1
        return rem
        
    def DistributeReads(self):
        confs = 0
        for rd in self.Reads.values():
            a_len = len(rd.alignments)
            rln = rd.length
            if a_len > 1 and not rd.matched:
                ks = list(rd.alignments.keys())
                als = list(rd.alignments.values())
                for i in range(a_len-1):
                    for j in range(i+1,a_len):
                        sc1 = self.Scaffolds[ks[i]]
                        sc2 = self.Scaffolds[ks[j]]
                        if sc1.CanRecieveNano(als[i], rln, sc2):
                            if sc2.CanRecieveNano(als[j], rln, sc1):
                                #print("Both Recievable!")
                                res = self.AttemptToReconcile(sc1,sc2,
                                                        als[i],als[j],rln)
                                if res:
                                    #print("Reconciled!")
                                    rd.matched += 1
                                    confs += 1
                
        return confs
                                
                            
                            
    def AttemptToReconcile(self, scaff1, scaff2, aln1, aln2, rlen):      
        def ExtensionLeg(aln, scj, rlen):
            extst = exted = 0
            
            if not aln1.rev and sc1j == 'R':
                extst = aln1.qend
                exted = rlen
            elif aln1.rev and sc1j == 'R':
                exted = aln1.qstart
                extst = 1
            elif not aln1.rev and sc1j == 'L':
                extst = 1
                exted = aln1.qstart
            elif aln1.rev and sc1j == 'L':
                extst = aln1.qend
                exted = rlen
                
            return [extst, exted]
        
        def Noverlap(aln1, aln2):
            
            if aln2.qstart <= aln1.qstart < aln2.qend:
                return False
            if aln2.qstart < aln1.qend <= aln2.qend:
                return False
            
            return True
        
        def Contains(ext, aln):
            if ext[0] <= aln.qstart < ext[1] and ext[0] < aln.qend <= ext[1]:
                return True
            return False
        
        def ResolveB(mid1, mid2, aln):
            scj = 'B'
            if mid1 < mid2:
                if aln.rev:
                    scj = 'L'
                else:
                    scj = 'R'
            else:
                if aln.rev:
                    scj = 'R'
                else:
                    scj = 'L'
                    
            return scj
            
        sensePass = False
        rangePass = False
        sizePass = False
        
        sc1j = scaff1.GetJoinStatus(scaff2.name)
        sc2j = scaff2.GetJoinStatus(scaff1.name)
        
        if not aln1.rev and not aln2.rev:
            if sc1j != sc2j or sc1j == 'B':
                sensePass = True
                
        elif aln1.rev != aln2.rev:
            if sc1j == sc2j or sc1j == 'B' or sc2j == 'B':
                sensePass = True
                
        elif aln1.rev and aln1.rev:
            if sc1j == sc2j == 'B':
                sensePass = True
        
        if(sensePass):
            ext1 = []
            ext2 = []
            
            cn1 = False
            cn2 = False
            
            if sc1j != 'B':
                ext1 = ExtensionLeg(aln1, sc1j, rlen)
                cn1 = Contains(ext1, aln2)
            else:
                cn1 = True
                
            if sc2j != 'B':
                ext2 = ExtensionLeg(aln2, sc2j, rlen)
                cn2 = Contains(ext2, aln1)
            else:
                cn2 = True
                           
            if cn1 and cn2 and Noverlap(aln1, aln2):
                rangePass = True
                
            if(rangePass):
                
                inds = sorted([aln1.qstart, aln1.qend, aln2.qstart, aln2.qend])
                mid1 = aln1.qstart + ((aln1.qend-aln1.qstart)/2)
                mid2 = aln2.qstart + ((aln2.qend-aln2.qstart)/2)
                
                qsep = abs(inds[1] - inds[2])
                
                clip = 0
                
                new_j1 = sc1j
                new_j2 = sc2j
                
                if sc1j == 'B':
                    new_j1 = ResolveB(mid1, mid2, aln1)
                if sc2j == 'B':
                    new_j2 = ResolveB(mid2, mid1, aln2)
                    
                sc1_clip = scaff1.GetTail(aln1, new_j1)
                sc2_clip = scaff2.GetTail(aln2, new_j2)
                
                bth_clip = sc1_clip + sc2_clip
                    
                clip = (bth_clip) - qsep
                
                gap = clip * -1
                tr1 = tr2 = 0
                if clip < TRIM_LIMIT:
                    sizePass = True
                    if clip > 0:
                        tr2 = int((sc2_clip / bth_clip) * clip)
                        tr1 = clip - tr2
                        gap = 0 
                    
                if(sizePass):
                    scaff1.AnchorNano(scaff2.name, sc1j, new_j1, tr1, gap, aln1.rev)
                    scaff2.AnchorNano(scaff1.name, sc2j, new_j2, tr2, gap, aln2.rev)
                    return 1
                    
        return 0
    
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
                    

class Agent:
    def __init__(self, scdict):
        self.Scaffolds = scdict
                
        
def HelpMsg():
    print ("\nNanoChrome - Confirming 10x with nanopore\n")
    print ("Input:\n")
    print ("ARG1: <genome.fa> A Genome file")
    print ("ARG2: <nc_table.tsv> output from chrome_candidates.py")
    print ("ARG3: <alns.sam> Mapped Longread library")
    print ("ARG4: <prefix> Unique prefix for outputs\n")
    print ("Output:\n")
    print ("<prefix>_nc_scaffolded.fa: Final Scaffolded Genome")
    print ("<prefix>_network_final.tsv: A network file for visualisation")
    print ("<prefix>_nodes_final.tsv: A node description file for visualisation\n\n")

def main(argv):

    # make sure there are at least three arguments
    if len(argv) >= 4:
        #try:
        the_graph = Graph()
        print("Scanning the genome...\n")
        the_graph.ScanGenome(argv[0])
        print("Reading the edge graph...\n")
        the_graph.ReadTable(argv[1])
        print("Pre-filtering edge graph...\n")
        the_graph.PreFilterEdges()
        the_graph.PreConfirmEdges()
        print("Loading long-read .paf...\n")
        the_graph.ReadPAF(argv[2])
        print("Filtering network...\n")
        the_graph.StripUnconfirmed()
        print("Agent navigating graph...\n")
        the_agent = Agent(the_graph.Scaffolds)
        
        print("Fixing in optimal scaffold paths\n")
        the_graph.StoreGenome(argv[0])
        print("Writing Final Assembly\n")
            
            
        #except:
        #    print("Error: ",sys.exc_info()[0]," <- this happened.")
        #finally:
        #   HelpMsg() 
        
    else:
        HelpMsg()
        print("You need to specify the four positional arguments\n")
        sys.exit(2)
 
 
if __name__ == '__main__':
    main(sys.argv[1:])