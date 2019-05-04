# -*- coding: utf-8 -*-

import sys
import pandas as pd
from collections import defaultdict

QUANTILE_THRESH = 0.999
EDGE_RATIO = 5
CONF_WEIGHT = 20
DEGREE_FILTER = 5
QUAL = 20
AUTO_GAP = 2000
TRIM_LIMIT = 500
FRAG_LEN = 50000
STRICT = 0

def ChooseDest(dct):
    def lstreset(jn):
        lst = []
        lst.append(jn)
        return lst
    
    if len(dct) > 0:
        best = next(iter(dct.values()))
        lst = lstreset(best)
        
        for ks,jn in dct.items():
            if jn.nano_confs > best.nano_confs:
                best = jn
                lst = lstreset(best)
            elif jn.nano_confs == best.nano_confs:
                if jn.auto_confs > best.auto_confs:
                    best = jn
                    lst = lstreset(best)
                elif jn.auto_confs == best.auto_confs:
                    if jn.weight > best.weight:
                        best = jn
                        lst = lstreset(best)
                    elif jn.nano_confs and jn.gap < best.gap:
                        best = jn
                        lst = lstreset(best)
                    elif jn.weight == best.weight:
                        lst.append(jn)
                            
        return lst
    else:
        return []

class JoinData:
    def __init__(self, wt, auto, nano):
        self.wt = wt
        self.auto = auto
        self.nano = nano

class Join:
    def __init__(self, other, weight, meanpos):
        self.other = other
        self.weight = weight
        self.nano_confs = 0
        self.auto_confs = 0
        self.orientation = 'F'
        self.gap = AUTO_GAP
        self.trim_amnt = 0
        self.mean_pos = meanpos
        
    def Reflect(self, scnm):       
        res, ref = self.other.HasFixedMatch(scnm)
        
        if res != 'N':
            return [res, ref]
        else:
            return [res, None]
    
    def __lt__(self,other):
        return self.weight < other.weight

class MetaScaffold:
    def __init__(self, scaffolds, joins, gaps, trims, jdat, entryd, exitd):
        self.Scaffolds = scaffolds
        self.JoinTypes = joins
        self.Gaps = gaps
        self.Trims = trims
        self.Jdata = jdat
        self.Entry = entryd
        self.Exit = exitd
        self.subsumed = False
        scln = sum([sc.length for sc in self.Scaffolds])
        
        self.length = scln + sum(self.Gaps) - sum(self.Trims)
        
        i = 0
        for sc in self.Scaffolds:
            if i == 0 or i == (len(self.Scaffolds)-1):
                sc.is_meta_terminal = True
            sc.SetInMeta(self)
            i += 1
            
    def GetSide(self, scaff):
        
        if scaff.name == self.Scaffolds[0].name:
            return 'L'
        elif scaff.name == self.Scaffolds[-1].name:
            return 'R'
        else:
            return ''
        
    def Mergelists(self, scaffs, exit_sc, gaps, trims, jd, joins, side, rev):
        if rev:
            scaffs = scaffs[::-1]
            gaps = gaps[::-1]
            trims = trims[::-1]
            jd = jd[::-1]
            joins = joins[::-1]
            
        if side == 'L':
            if scaffs:
                self.Scaffolds = scaffs + self.Scaffolds
            self.Entry = exit_sc
            self.Gaps = gaps + self.Gaps
            self.Trims = trims + self.Trims
            self.Jdata = jd + self.Jdata
            self.JoinTypes = joins + self.JoinTypes
        
        if side == 'R':
            if scaffs:
                self.Scaffolds = self.Scaffolds + scaffs
            self.Exit = exit_sc
            self.Gaps = self.Gaps + gaps
            self.Trims = self.Trims + trims
            self.Jdata = self.Jdata + jd
            self.JoinTypes = self.JoinTypes + joins
        
    def Subsume(self, term_sc, meta, jn_side, my_side):
        meta.subsumed = True
        #First connect the two edge scaffolds
        
        sc1 = self.Scaffolds[0] if my_side == 'L' else self.Scaffolds[-1]
        sc2 = meta.Scaffolds[0] if jn_side == 'L' else meta.Scaffolds[-1]
        
        ch1, jn1 = sc1.HasFixedMatch(sc2.name)
        ch2, jn2 = sc2.HasFixedMatch(sc1.name)
        
        if ch1 != 'N' and ch2 != 'N':
            gaps = [jn1.gap]
            join_ty = [jn1.orientation, jn2.orientation]
            trims = [jn1.trim_amnt, jn2.trim_amnt]
            join_dat = [JoinData(jn1.weight, jn1.auto_confs, jn1.nano_confs)]
            exdm = 'P'
            
            rev_mini = True if my_side == 'L' else False
            self.Mergelists([], exdm, gaps, trims,
                            join_dat, join_ty, my_side, rev_mini)
            
            #Then merge two metascaffs here (mostly their att-lists)
            rev = True if jn_side == my_side else False
            ex_sc = meta.Exit if jn_side == 'L' else meta.Entry
            self.Mergelists(meta.Scaffolds, ex_sc, meta.Gaps, meta.Trims,
                            meta.Jdata, meta.JoinTypes, my_side, rev)
            return True
#        else:
#            print ("Not finding matches despite traversal!\n")
        
        return False
            
    def Append(self, prex):
        subbed = 0
        for sc in prex.pre_meta.scaffs:
            sc.in_extend = True
            sc.in_meta = True
        
        side = prex.side
        pm = prex.pre_meta
        
        if len(pm.scaffs) > 1:
            rev = True if side == 'L' else False
            self.Mergelists(pm.scaffs[1:], pm.exit, pm.gaps, pm.trims,
                            pm.jdat, pm.joins, side, rev)
            
        if prex.pre_meta.meta_terminal:
            edge_to_join = prex.pre_meta.scaff_last
            met_to_join = edge_to_join.meta_ref
            join_side = met_to_join.GetSide(edge_to_join)
            if not met_to_join.subsumed:
                if self.Subsume(edge_to_join, met_to_join, join_side, side):
                    subbed = 1
            
        scln = sum([sc.length for sc in self.Scaffolds])
        
        self.length = scln + sum(self.Gaps) - sum(self.Trims)
        
        for s in self.Scaffolds:
            s.is_meta_terminal = False
            
        self.Scaffolds[0].is_meta_terminal = True
        self.Scaffolds[-1].is_meta_terminal = True
        
        self.StripTerminalInners()
        
        return subbed
        
    
    def StripTerminalInners(self):
        def strin(sc, ch):
            if ch == 'R':
                sc.QuarantineIt(L=True, R=False, B=True)
                
            if ch == 'L':
                sc.QuarantineIt(L=False, R=True, B=True)
            
        sc1 = self.Scaffolds[0]
        sc2 = self.Scaffolds[-1]
        
        strin(sc1, self.Entry)
        strin(sc2, self.Exit)
        
            
    def Name(self, name):
        self.Name = name
    
    def PrintNet(self, nodes, net):
        for sc in self.Scaffolds:
            nodes.write(sc.name + "_L\t" + str(sc.length) + "\n")
            nodes.write(sc.name + "_R\t" + str(sc.length) + "\n")
            net.write(sc.name + "_L\t" + sc.name + "_R\t2000\t3\n")
            
        for i in range(len(self.Jdata)):
            sc1nm = self.Scaffolds[i].name + "_R"
            sc2nm = self.Scaffolds[i+1].name + "_L"
            stub = sc1nm + "\t" + sc2nm + "\t" + str(self.Jdata[i].wt)
            auto = self.Jdata[i].auto
            nano = self.Jdata[i].nano
            if auto:
                net.write(stub + "\t1\n")
            if nano:
                net.write(stub + "\t2\n")
            
    
    def PrintIt(self, handle):
        old_chars = "ACGT"
        replace_chars = "TGCA"
        revcomp = str.maketrans(old_chars,replace_chars)
        
        seqs = []
        main_seq = ""
        for sc in self.Scaffolds:
            seqs.append(sc.seq)
        
        lastrev = False
        
        if len(seqs) > 0:
            handle.write(">" + self.Name + "\n")
            for i in range(len(seqs)):
                tr1 = tr2 = 0
                j1 = 'F'
                if i > 0:
                    tr1 = self.Trims[(i*2)-1]
                    j1 = self.JoinTypes[(i*2)-1]
                if i < (len(seqs) - 1):
                    tr2 = self.Trims[i*2]
                tr2 = len(seqs[i]) - tr2
                
                flip = False
                
                if lastrev and j1 == 'F':
                    flip = True
                elif not lastrev and j1 == 'R':
                    flip = True
                
                if flip:
                    seqs[i] = seqs[i].translate(revcomp)[::-1]
                
                main_seq += seqs[i][tr1:tr2]
                
                if i < (len(seqs) - 1):
                    main_seq += 'N' * self.Gaps[i]
                    lastrev = flip
                    
        mln = len(main_seq)
        inc = 0
        while(inc < (mln - 100)):
            subs = main_seq[inc:(inc + 100)]
            handle.write(subs + "\n")
            inc += 100
        
        rem = mln % 100
        if rem > 0:
            subs = main_seq[(mln - rem):mln]
            handle.write(subs + "\n")

class Scaffold:
    def __init__(self, name, length):
        self.name = name
        self.length = length
        self.in_meta = False
        self.in_loop = False
        self.in_extend = False
        self.seq = ""
        self.is_meta_terminal = False
        self.meta_ref = None
        
        self.Left = {}
        self.Right = {}
        self.Both = {}
        
        self.Quarantine = defaultdict(list)
        
    def ShuntBoths(self):
        if self.Both:
            if self.Left and not self.Right:
                self.DegreeFilter(l=False,r=False,b=True)
                self.Right = self.Both
                self.Both = {}
            if not self.Left and self.Right:
                self.DegreeFilter(l=False,r=False,b=True)
                self.Left = self.Both
                self.Both = {}
    
    def QuarantineIt(self, L=True, R=True, B=True):
        def Qdict(dct, ch):
            for jn,vl in dct.items():
                self.Quarantine[ch].append(vl)
                vl.other.QuarantineOther(self.name)
            
            return {}
        
        if L:
            self.Left = Qdict(self.Left, 'L')
        
        if R:
            self.Right = Qdict(self.Right, 'R')
        
        if B:
            self.Both = Qdict(self.Both, 'B')
    
    def DeQuarantine(self):
        def DeQr(dct, ch):
            if ch in self.Quarantine:
                for jn in self.Quarantine[ch]:
                    if not jn.other.in_meta or jn.other.is_meta_terminal:
                        dct[jn.other.name] = jn
        
        if not self.in_meta or self.is_meta_terminal:
            DeQr(self.Left, 'L')
            DeQr(self.Right, 'R')
            DeQr(self.Both, 'B')
            
    def SetInMeta(self, met):
        self.QuarantineIt()            
        self.in_meta = True
        self.meta_ref = met
        
    def PrintIt(self, handle, newid):
        mln = len(self.seq)
        inc = 0
        handle.write(">original_seq_" + str(newid) + "\n")
        while(inc < (mln - 100)):
            subs = self.seq[inc:(inc + 100)]
            handle.write(subs + "\n")
            inc += 100
        
        rem = mln % 100
        if rem > 0:
            subs = self.seq[(mln - rem):mln]
            handle.write(subs + "\n")
            
    def IsBest(self, jn, nxt):
        def ReflectDest(dct, b1):
            if len(dct) > 0:
                best = b1
                result = 2
                for ks,jn in dct.items():
                    if jn.nano_confs > best.nano_confs:
                        best = jn
                        result = 0
                    elif jn.nano_confs == best.nano_confs:
                        if jn.auto_confs > best.auto_confs:
                            best = jn
                            result = 0
                        elif jn.auto_confs == best.auto_confs:
                            if jn.weight > best.weight:
                                best = jn
                                result = 0
                            elif jn.nano_confs and jn.gap < best.gap:
                                best = jn
                                result = 0
                            elif jn.weight == best.weight:
                                result = 1
                return [best, result]
            else:
                return [None, 0]
            
        res = 0
        
        if nxt == 'L':
            if len(self.Left) == 1:
                res = 2
            else:
                jn2, res = ReflectDest(self.Left, jn)
                
        if nxt == 'R':
            if len(self.Right) == 1:
                res = 2
            else:
                jn2, res = ReflectDest(self.Right, jn)
            
        return res
    
        
    def IsEntryPoint(self):
        def HasRBJ(dct):
            lst = ChooseDest(dct)
            for jn in lst:
                next_entry, jn2 = jn.Reflect(self.name)
                rbj_result = jn.other.IsBest(jn2, next_entry) if jn2 else 0
                if rbj_result:
                    return True
            else:
                return False
        left_rbj = HasRBJ(self.Left)
        right_rbj = HasRBJ(self.Right)
        
        if left_rbj and not right_rbj:
            return 'L'
        if right_rbj and not left_rbj:
            return 'R'
        
        return ''
    
        
    def RemoveTakenJoins(self):
        def stmeta(dct, ch):
            ks = list(dct.keys())
            for k in ks:
                if dct[k].other.in_meta:
                    self.Quarantine[ch].append(dct[k])
                    del dct[k]
                    
        if not self.in_meta:
            stmeta(self.Left, 'L')
            stmeta(self.Right, 'R')
            #stmeta(self.Both, 'B')
    
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
            jn.nano_confs += 1
            if rev:
                jn.orientation = 'R'
                
        self.SwitchDict(old_mode, new_mode, othnm)
        
        if new_mode == 'L':
            UpdateSwitch(self.Left[othnm], trim, gap, rev)
            
        elif new_mode == 'R':
            UpdateSwitch(self.Right[othnm], trim, gap, rev)
        
    
    def HasFixedMatch(self, scfnm):
        if scfnm in self.Left:
            return ['L', self.Left[scfnm]]
        elif scfnm in self.Right:
            return ['R', self.Right[scfnm]]
        else:
            return ['N', None]
        
             
    def RunAutoConf(self, fraglen):
        def CheckEntry(tail1, vl, ch1, ch2):
            res, v2 = vl.other.HasFixedMatch(self.name)
            if res == 'L':
                tail2 = v2.mean_pos
                remains = int(fraglen - tail1 - tail2)
                vl.orientation = ch1
                vl.auto_confs = 1
                vl.gap = max([TRIM_LIMIT, remains])
                return 1
            if res == 'R':
                tail2 = vl.other.length - v2.mean_pos
                remains = int(fraglen - tail1 - tail2)
                vl.orientation = ch2
                vl.auto_confs = 1
                vl.gap = max([TRIM_LIMIT, remains])
                return 1
            return 0
        
        cnf = 0
        allc = 0
        for ks,vl in self.Left.items():
            if not vl.auto_confs and not vl.nano_confs:
                allc += 1
                tail1 = vl.mean_pos
                cnf += CheckEntry(tail1, vl, 'R', 'F')
            
        for ks,vl in self.Right.items():
            if not vl.auto_confs and not vl.nano_confs:
                allc += 1
                tail1 = self.length - vl.mean_pos
                cnf += CheckEntry(tail1, vl, 'F', 'R')
        
        allc += len(self.Both)
        
        return [allc, cnf]
        
    def AddJoin(self, other, weight, mpos, mode):
        jn = Join(other, weight, mpos)
        
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
            
    def DegreeFilter(self, l=True,r=True,b=False):
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
                hval = ln - DEGREE_FILTER
                k, c = SplitDict(dct)
                replacement = k
                FilterOut(c)
                    
            return replacement, hval
        
        if l:
            ret, hv = FilterDict(self.Left)
            self.Left = ret
            hval += hv
        
        if r:
            ret, hv = FilterDict(self.Right)
            self.Right = ret
            hval += hv
        
        if b:
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
            
    def QuarantineOther(self, name):
        if name in self.Right:
            self.Quarantine['R'].append(self.Right[name])
            del self.Right[name]
        
        if name in self.Left:
            self.Quarantine['L'].append(self.Left[name])
            del self.Left[name]
            
        if name in self.Both:
            self.Quarantine['B'].append(self.Both[name])
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
    def __init__(self, flen):
        self.Scaffolds = {}
        self.Reads = {}
        self.PreN50 = 0
        self.PostN50 = 0
        self.FragLen = flen
        
    def QuarantineEverything(self):
        for sc,vl in self.Scaffolds.items():
            vl.QuarantineIt()
        
        for sc,vl in self.Scaffolds.items():
            vl.DeQuarantine()
        
    def ResetScaffs(self):
        for sc,vl in self.Scaffolds.items():
            vl.in_meta = False
            vl.in_loop = False
        
    def GetSeq(self, othnm):
        if othnm in self.Scaffolds:
            return self.Scaffolds[othnm].seq
        
    def RunN50(self, lns):
        lns.sort()
        
        sz = sum(lns)/2
        accu = 0
        n50 = 0
        oldL = 0
        for l in lns:
            if accu + l > sz:
                a = (sz - accu) / l
                b = l - oldL
                n50 = oldL + (b * a)
                break;
            else:
                accu += l
                oldL = l
        
        return n50
    
    def FinalStats(self, metas):
        lns = []
        
        for ks,vl in self.Scaffolds.items():
            if not vl.in_meta:
                lns.append(vl.length)
        
        for m in metas:
            lns.append(m.length)
            
        ln = len(lns)
        sm = sum(lns)
            
        return [self.RunN50(lns), ln, sm]
    
    def CalcN50(self, dct):
        
        lns = [0] * len(dct)
        i = 0
        for ks,vl in dct.items():
            lns[i] = vl.length
            i += 1
        
        n50 = self.RunN50(lns)
        
        return n50
        
        
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
            allc1, cnf1 = scs.RunAutoConf(self.FragLen)
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
                    nm = nm.split(' ')[0]
                    if len(oldID) > 0 and len(seq) > 15:
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
        
        self.PreN50 = self.CalcN50(self.Scaffolds)
        
        print("Read Genome with " + str(glen) + " bases")
        print("Found " + str(snum) + " separate contigs")
        print("Assembly N50 is: %.2f\n" % self.PreN50)
        
    def ReadTable(self, tnam):
        def AddRow(row, members):
            sc1 = self.Scaffolds[row.CNTG1]
            sc2 = self.Scaffolds[row.CNTG2]
            members[row.CNTG1] = 1
            members[row.CNTG2] = 1
            sc1.AddJoin(sc2, row.WEIGHT, row.MEAN1, row.TYPE1)
            sc2.AddJoin(sc1, row.WEIGHT, row.MEAN2, row.TYPE2)
            
        def GetInc(t1, t2, tst):
            rv = 0
            if t1 == tst:
                rv += 1
                
            if t2 == tst:
                rv += 1
                
            return rv
        
        jinc = linc = rinc = binc = 0
        
        tbl = pd.read_csv(tnam, sep="\t")
        tbl['CNTG1'] = tbl['CNTG1'].astype(str)
        tbl['CNTG2'] = tbl['CNTG2'].astype(str)
        
        print("Actual join count: " + str(tbl.shape[0]))
        global CONF_WEIGHT
        global QUANTILE_THRESH
        
        sln = len(self.Scaffolds)
        act_ratio = 0
        
        print("Adjusting quantile filter...\n")
        
        tblA = tbl[tbl['TYPE1'] != 'B']
        tblA = tblA[tblA['TYPE2'] != 'B']
        
        tblB = tbl[tbl['TYPE1'] == 'B']
        tblB = tblB[tblB['TYPE2'] == 'B']
        
        tblX = tbl[tbl['TYPE1'] != 'B']
        tblX = tblX[tblX['TYPE2'] == 'B']
        
        tblY = tbl[tbl['TYPE1'] == 'B']
        tblY = tblY[tblY['TYPE2'] != 'B']
        
        tbl = tblX.append(tblY)
        
        while QUANTILE_THRESH > 0.7 and act_ratio < EDGE_RATIO and CONF_WEIGHT >= 5: 
            CONF_WEIGHT = tbl.WEIGHT.quantile(QUANTILE_THRESH)
            print("%.2f quantile filter; weight threshold: %.2f"
                  % (QUANTILE_THRESH, CONF_WEIGHT))
            tbl2 = tbl.loc[tbl.WEIGHT > CONF_WEIGHT]
            ln2 = tbl2.shape[0]
            act_ratio = ln2/sln
            QUANTILE_THRESH -= 0.005
            
        tbl2 = tbl.loc[tbl.WEIGHT > CONF_WEIGHT]
        tblA = tblA.loc[tblA.WEIGHT > CONF_WEIGHT]
        print("\nPost filter partials to load: " + str(tbl2.shape[0]))
        print("Chromium direct tangle to load: " + str(tblA.shape[0]))
        print("Small Frag Nano-bait to load: " + str(tblB.shape[0]) + "\n")
        
        membership = {}
        
        CONF_WEIGHT = int(CONF_WEIGHT)
        
        for ind,row in tbl2.iterrows():      
            jinc += 2
            rinc += GetInc(row.TYPE1,row.TYPE2, 'R')
            linc += GetInc(row.TYPE1,row.TYPE2, 'L')
            binc += GetInc(row.TYPE1,row.TYPE2, 'B')
            
            AddRow(row, membership)
            
        for ind,row in tblB.iterrows():         
            jinc += 2
            binc += GetInc(row.TYPE1,row.TYPE2, 'B')
            AddRow(row, membership)
            
        for ind,row in tblA.iterrows():         
            jinc += 2
            rinc += GetInc(row.TYPE1,row.TYPE2, 'R')
            linc += GetInc(row.TYPE1,row.TYPE2, 'L')
            AddRow(row, membership)
            
        lnmb = len(membership)
            
        print ("Loaded " + str(jinc) + " partial joins")
        print ("With " + str(linc) + " left arks")
        print ("With " + str(rinc) + " right arks")
        print ("With " + str(binc) + " ambiguous arks")
        print ("Shared by " + str(lnmb) + " scaffolds\n")
        
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
                                res = self.AttemptToReconcile(sc1,sc2,
                                                        als[i],als[j],rln)
                                if res:
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
                    nm = nm.split(' ')[0]
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
        self.Scaffolds = list(scdict.values())
        self.Metas = []
        
    class PreExtend:
        def __init__(self, prem, met, side):
            self.side = side
            self.pre_meta = prem
            self.meta_scaff = met
        
        def CheckExtend(self):
            if len(self.pre_meta.scaffs) > 1:
                for sc in self.pre_meta.scaffs:
                    if sc.in_extend:
                        return False
                else:
                    return True
            elif self.pre_meta.meta_terminal:
                if not self.pre_meta.scaff_last.meta_ref.subsumed:
                    return True
            return False
        
        def __lt__(self, other):
            return self.pre_meta.__lt__(other.pre_meta)
        
    class PreMeta:
        def __init__(self, entry):
            self.scaffs = []
            self.joins = []
            self.gaps = []
            self.trims = []
            self.jdat = []
            self.entry = entry
            self.exit = entry
            self.loop_terminal = False
            self.meta_terminal = False
            self.scaff_last = None
            self.last_entry = ''
            
        def LoopCheck(self):
            for sc in self.scaffs:
                if sc.in_loop:
                    return False
                
            return True
        
        def Shatter(self, rbj_thresh):
            def NonRBHQuarantine(dct, sc, ch, rbj):
                ndict = {}
                for jn,vl in dct.items():
                    nxt, jn2 = vl.Reflect(sc.name)
                    rbj_result = vl.other.IsBest(jn2, nxt)
                    if rbj_result < rbj:
                        sc.Quarantine[ch].append(vl)
                        vl.other.QuarantineOther(sc.name)
                    else:
                        ndict[jn] = vl
                        
                return ndict
                        
            def CheckQr(ch, sc):
                dif = 0
                if ch == 'L':
                    lns1 = len(sc.Left)
                    sc.Left = NonRBHQuarantine(sc.Left, sc, ch, rbj_thresh)
                    dif = lns1 - len(sc.Left)
                if ch == 'R':
                    lns1 = len(sc.Right)
                    sc.Right = NonRBHQuarantine(sc.Right, sc, ch, rbj_thresh)
                    dif = lns1 - len(sc.Right)
                    
                return dif
            
            diff = 0
            for sc in self.scaffs:
                sc.in_loop = True
                
            sc1 = self.scaffs[0]
            sc2 = self.scaffs[-1]
            
            diff += CheckQr(self.entry, sc1)
            diff += CheckQr(self.exit, sc2)
            
            return diff
              
            
        def AttemptClosure(self):
            for s in self.scaffs:
                if s.in_meta:
                    return None
            
            met = MetaScaffold(self.scaffs, self.joins,
                               self.gaps, self.trims, self.jdat,
                               self.entry, self.exit)
            
            return met
            
        def __lt__(self, other):
            ln = len(self.scaffs)
            ln2 = len(other.scaffs)
            
            if ln != ln2:
                return ln < ln2
            else:
                ln = sum([sc.length for sc in self.scaffs])
                ln2 = sum([sc.length for sc in other.scaffs])

                return ln < ln2
            
    def ResolveNewlyDirectionalJoins(self):
        for sc in self.Scaffolds:
            sc.ShuntBoths()
        
    def StripGraph(self):
        sc2 = []
        for sc in self.Scaffolds:
            sc.RemoveTakenJoins()
            if not sc.in_meta:
                sc2.append(sc)
                
        self.Scaffolds = sc2
        
    def StripMetas(self):
        mt2 = []
        for m in self.Metas:
            if not m.subsumed:
                mt2.append(m)
                
        self.Metas = mt2
        
    def Traverse(self, scaff, entry):
        entry_point = entry
        
        journey = {}
        pm = self.PreMeta(entry_point)
        
        can_move = True
        
        while(can_move):
            journey[scaff.name] = scaff
            pm.scaffs.append(scaff)
            
            jn = None
            lstx = []
            if entry_point == 'L':
                lstx = ChooseDest(scaff.Right)
            if entry_point == 'R':
                lstx = ChooseDest(scaff.Left)
                
            if not lstx:
                can_move = False
            else:
                rbj_result = rbjo = 0
                next_entry = nxto = ''
                jn2 = jn2o = jno = None
                
                for jnx in lstx:
                    jn = jnx
                    next_entry, jn2 = jn.Reflect(scaff.name)
                    rbj_result = jn.other.IsBest(jn2, next_entry) if jn2 else 0
                    
                    if rbj_result == 2:
                        break
                    
                    if rbj_result == 1:
                        jn2o = jn2
                        jno = jn
                        nxto = next_entry
                        rbjo = 1
                        
                    if STRICT and rbj_result != 2:
                        rbj_result = 0
                        break
                        
                strictpass = True if not STRICT else (len(lstx) < 2) 
                
                if not rbj_result and rbjo and strictpass:
                    jn2 = jn2o
                    jn = jno
                    rbj_result = rbjo
                    next_entry = nxto
                
                inj = jn.other.name in journey
                
                if inj:
                    pm.loop_terminal = True
                elif jn.other.in_meta:
                    pm.meta_terminal = True
                    pm.scaff_last = jn.other
                    pm.last_entry = next_entry
                    
                if not jn2 or jn.other.in_meta or not rbj_result or inj or jn.other.in_extend:
                    can_move = False                    
                else:
                    scaff = jn.other
                    entry_point = next_entry
                    pm.exit = 'L' if entry_point == 'R' else 'R'
                    pm.joins.append(jn.orientation)
                    pm.joins.append(jn2.orientation)
                    pm.gaps.append(jn.gap)
                    pm.trims.append(jn.trim_amnt)
                    pm.trims.append(jn2.trim_amnt)
                    pm.jdat.append(JoinData(jn.weight, 
                                            jn.auto_confs, jn.nano_confs))
                
        return pm
    
    def LoopShatter(self, rbj_thresh, limit):
        sh_cnt = 0
        ed_cnt = 0
        sc_cnt = 0
        pms = []
        
        for sc in self.Scaffolds:
            if sc.Left and sc.Right and not sc.in_meta:
                pma = self.Traverse(sc, 'L')
                pmb = self.Traverse(sc, 'R')
                
                if pma.loop_terminal and pmb.loop_terminal:
                    if len(pma.scaffs) > 1:
                        pms.append(pma)
                    if len(pmb.scaffs) > 1:
                        pms.append(pmb)
                    
        pms = sorted(pms, reverse=True)
        for p in pms:
            if p.LoopCheck():
                res = p.Shatter(rbj_thresh)
                ed_cnt += res
                sc_cnt += 1 if res else 0
                sh_cnt += 1
                if sc_cnt >= limit:
                    break
                
        return [sh_cnt, ed_cnt, sc_cnt]
                
        
    def RunPaths(self):
        def RunDisentangler(rbj, iteration_limit):
            for sc in self.Scaffolds:
                sc.in_loop = False
                
            print("Running disentangler..")
            
            shat, edg, succ = self.LoopShatter(rbj, iteration_limit)
            print("Ran on " + str(shat) + " tangles, rbj: " + str(rbj))
            print("..picked open " + str(succ) + " loops")
            print("..quarantining " + str(edg) + " joins\n")
            return [succ, shat]
            
        metas_built = 1
        metas_kept = 1
        total_metas = 0
        inc = 0
        prems = []
        lsht = 0
        
        tscs = len(self.Scaffolds)
        inmetas = 0
        succ = 0
        shat = 0
        
        iteration_limit = int(tscs/250)
        
        while(metas_kept > 0):
            metas_built = 0
            metas_kept = 0
            succ = 0
            shat = 0
            prems = []
            inc += 1
            self.StripGraph()
            
            if inc > 0:
                rbj =  int(total_metas/(tscs/25))
                succ, shat = RunDisentangler(rbj, iteration_limit)
                lsht += succ
            
            for sc in self.Scaffolds:
                res = sc.IsEntryPoint()
                if res:
                    res = 'L' if res == 'R' else 'R'
                    npm = self.Traverse(sc, res)
                    
                    if len(npm.scaffs) > 1:
                        prems.append(npm)
                        metas_built += 1
            
            if len(prems) > 1:
                prems = sorted(prems, reverse=True)
            
            for pm in prems:
                met = pm.AttemptClosure()
                if met:
                    metas_kept += 1
                    total_metas += 1
                    inmetas += len(pm.scaffs)
                    met.Name("Nanochrome_scaffold_" + str(total_metas))
                    self.Metas.append(met)
                    if metas_kept > iteration_limit:
                        break    
            
            print("Meta-scaffolding iteration: " + str(inc))
            print("..new scaffolds created: " + str(metas_built))
            print("..of which kept: " + str(metas_kept) + "\n")
            
            if not metas_kept and succ:
                metas_kept = succ
            elif not metas_kept and not succ and shat and rbj < 2:
                if inc > 0:
                    succ, shat = RunDisentangler(4, iteration_limit)
                    lsht += succ
                    
            if not metas_kept and succ:
                metas_kept = succ
                
        print ("Grand Total of: " + str(total_metas) + " metascaffolds built")
        print ("Unpicked " + str(lsht) + " tangles in graph")
        print ("Of " + str(tscs) + " scaffolds in genome, " 
               + str(inmetas) + " were used\n")
        
    def ExtendMetas(self):
        for m in self.Metas:
            m.StripTerminalInners()
        
        total_extend = 0
        scaffs_included = 0
        metas_joined = 0
        metas_extend = 1
        extend_built = 0
        lsht = 0
        inc = 0
        
        while(metas_extend):
            metas_extend = 0
            extend_built = 0
            extends = []
            inc += 1
            
            self.StripMetas()
            
            for sc in self.Scaffolds:
                sc.in_loop = False
            
            if inc > 1:
                print("Running disentangler..")
                shat, edg, succ = self.LoopShatter(inc, 50)
                lsht += succ
                print("Ran on " + str(shat) + " tangles")
                print("..picked open " + str(succ) + " loops")
                print("..quarantining " + str(edg) + " joins\n")
                
            for m in self.Metas:
                scL = m.Scaffolds[0]
                scR = m.Scaffolds[-1]
                if scL and scR:
                    resL = scL.IsEntryPoint()
                    resR = scR.IsEntryPoint()
                    
                    if resL:
                        resL = 'L' if resL == 'R' else 'R'
                        pmL = self.Traverse(scL, resL)
                        
                        if len(pmL.scaffs) > 1:
                            prex = self.PreExtend(pmL, m, 'L')
                            extends.append(prex)
                            extend_built += 1
                            
                    if resR:
                        resR = 'L' if resR == 'R' else 'R'
                        pmR = self.Traverse(scR, resR)
                        
                        if len(pmR.scaffs) > 1:
                            prex = self.PreExtend(pmR, m, 'R')
                            extends.append(prex)
                            extend_built += 1
                        
            if len(extends) > 1:
                extends = sorted(extends, reverse=True)
            
            for ex in extends:
                if ex.CheckExtend():
                    metas_joined += ex.meta_scaff.Append(ex)
                    metas_extend += 1
                    total_extend += 1
                    scaffs_included += (len(ex.pre_meta.scaffs) - 1)
            
            print ("Meta-scaff extention iteration: " + str(inc))
            print ("..extensions made: " + str(extend_built))
            print ("..of which kept: " + str(metas_extend) + "\n")
            
            if not metas_extend and succ:
                metas_extend = 1
            
    
        fln = len(self.Metas)
        print ("Grand Total of: " + str(fln) + " metascaffolds remain")
        print ("To which " + str(total_extend) + " extensions were made")
        print ("Unpicked " + str(lsht) + " tangles in graph")   
        print ("Newly included contigs: " + str(scaffs_included))
        print ("Meta-scaffolds merged: " + str(metas_joined) + "\n")
        
        
    def PrintNetwork(self, pref):
        nodes = pref + "_nodes_final.tsv"
        net = pref + "_net_final.tsv"
        
        with open(nodes, "w") as ndfile:
            with open(net, "w") as ntfile:
                ndfile.write("NODE\tLENGTH\n")
                ntfile.write("SOURCE\tTARGET\tWEIGHT\tTYPE\n")
                for m in self.Metas:
                    m.PrintNet(ndfile, ntfile)
        
        print("Wrote " + nodes + " (network nodes)")
        print("Wrote " + net + " (network edges)\n")
                   
    def PrintMetas(self, pref, sc_resid):
        genm = pref + "_final_genome.fa"
        pm = ps = 0
        with open(genm, "w") as genomefile:
            for m in self.Metas:
                pm += 1
                m.PrintIt(genomefile)
            
            for ks,sc in sc_resid.items():
                if not sc.in_meta:
                    ps += 1
                    sc.PrintIt(genomefile, ps)
        print("Printed " + str(pm) + " metascaffs")
        print("Printed " + str(ps) + " scaffolds")
        print("Wrote " + genm + " successfully!\n")
       
def HelpMsg():
    print ("\nNanochrome - confirming 10x with nanopore\n")
    print ("Input:\n")
    print ("ARG1: <genome.fa> A Genome file")
    print ("ARG2: <nc_table.tsv> output from chrome_candidates.py")
    print ("ARG3: <alns.paf> Mapped Longread library")
    print ("ARG4: <integer> Fragment Length (x10 library)")
    print ("ARG5: <integer> Tangle leniency [5-15]")
    print ("ARG6: <0/1> Strict Mode (if re-run) [0]")
    print ("ARG7: <prefix> Unique prefix for outputs\n")
    print ("Output:\n")
    print ("<prefix>_nc_scaffolded.fa: Final Scaffolded Genome")
    print ("<prefix>_network_final.tsv: A network file for visualisation")
    print ("<prefix>_nodes_final.tsv: A node description file for visualisation\n\n")

def main(argv):

    # make sure there are seven arguments
    if len(argv) == 7:
        try:
            global FRAG_LEN
            global EDGE_RATIO
            global STRICT
            STRICT = int(argv[5])
            EDGE_RATIO = int(argv[4])
            FRAG_LEN = int(argv[3])
            the_graph = Graph(FRAG_LEN)
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
            the_agent = Agent(the_graph.Scaffolds)
            the_agent.ResolveNewlyDirectionalJoins()
            the_graph.PreConfirmEdges()
            the_graph.StripUnconfirmed()
            print("Agent navigating graph...\n")
            the_graph.ResetScaffs()
            the_agent.RunPaths()
            print("Agent re-navigating scaffold terminals...\n")
            the_graph.QuarantineEverything()
            the_agent.ExtendMetas()
            print("Fixing in optimal scaffold paths\n")
            the_agent.PrintNetwork(argv[6])
            n50, cnt, sz = the_graph.FinalStats(the_agent.Metas)
            
            print("The new genome size: " + str(sz))
            print("..with contig count: " + str(cnt))
            print("..and an n50 of: %.2f\n" % n50)
            
            print("Writing final assembly\n")
            the_graph.StoreGenome(argv[0])
            the_agent.PrintMetas(argv[6], the_graph.Scaffolds)
            
        except:
            print("Error: ",sys.exc_info()[0]," <- this happened.")
            HelpMsg()
    else:
        HelpMsg()
        print("You need to specify the four positional arguments\n")
        #sys.exit(2)
 
 
if __name__ == '__main__':
    main(sys.argv[1:])