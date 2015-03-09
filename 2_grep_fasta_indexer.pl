#!/usr/bin/env python2.7

import sys




if len(sys.argv) is not 3:
    print "Usage: <draft_vars.txt> <draft_vars.fasta>"
    sys.exit()
in_1 = open(sys.argv[1], "rU")
out_1 = open(sys.argv[2], "w")




count = 0 

with open(sys.argv[1]) as f:
  
    
    master = "".join(line.strip() for line in f)
    seq = master.split('--')
    for thing in seq:
        count = count +1
        print >>out_1,">" +str(count)
        print >>out_1,thing
    
    
    
    
    
    
    
out_1.close()