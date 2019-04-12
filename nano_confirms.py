# -*- coding: utf-8 -*-

import sys
import numpy as np

DEGREE_FILTER = 5

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
            
            
    def DegreeFilter(self):
        def FilterOut(dct):
            for jns in dct.values():
                jns.other.KillOther(self.name)
                
        if len(self.Left) > DEGREE_FILTER:
            FilterOut(self.Left)
            self.Left = {}
            
        if len(self.Right) > DEGREE_FILTER:
            FilterOut(self.Right)
            self.Right = {}
            
        if len(self.Both) > DEGREE_FILTER:
            FilterOut(self.Both)
            self.Both = {}
            
            
    def KillOther(self, name):
        if name in self.Right:
            del self.Right[name]
        
        if name in self.Left:
            del self.Left[name]
            
        if name in self.Both:
            del self.Both[name]
        

class Graph:
    def __init__(self):
        self.Scaffolds = {}
    
    def ReadTable(self, tnam):
        
        with open(tnam, "r") as tblin:
            next(tblin)
            for line in tblin:
                line.rstrip('\r\n')
                lns = line.split('\t')
                
                if lns[0] not in self.Scaffolds:
                    self.Scaffolds[lns[0]] = Scaffold(lns[0])

                
def HelpMsg():
    print ("NanoChrome - Confirming 10x with nanopore\n\n")
    print ("ARG1: <genome.fa> A Genome file\n")
    print ("ARG2: <alns.sam> Mapped Longread library\n")
    print ("ARG3: <nc_table.tsv> output from chrome_candidates.py")
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
            the_graph.ReadTable()
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