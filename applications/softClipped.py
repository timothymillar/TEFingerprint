import sys

"HW-ST994:243:C0T7UACXX:7:1110:17359:83679       99      CAULIV1-chr15_1 1       40      30S71M  =       68      168     GCTATTAGGCTCTAGATTAAAACCCTAAAATGAGTTCAATTTTCAGAACCATGTTTATGAGATCCACTGATTCTCAAGATACTATTGTAAACAGTGAAAAG   @?@FDEFDHDBAFHGGDGIEHIF@FCBGHDHIEEGHDFCDDGHD9?FGGHHGGCIG@F>GCDBF@F=B>=@C>C=.=3@EAEH@AEEE@;BCD>;ACE@>C   NM:i:0  MD:Z:71 AS:i:71 XS:i:73"


TYPES = ["S", "M", "D", "I"]
fractionation = []

O = open(sys.argv[2], 'w')

with open(sys.argv[1], 'r') as B:
#split up the line
    for line in B:
        line = line.split()
        CIGAR = line[5]
        sequence = line[9]
        quality = line[10]
        fraction = ""
#test if any softclipping occured 
        if "S" in CIGAR:
#chop by segment type
            for L in CIGAR:
                if L not in TYPES:
                    fraction += L
                else:
                    fraction += L
                    fractionation.append(fraction)
                    fraction = ""
            #Looks specifically at first or last elements
            firstEl = fractionation[0]
            lastEl = fractionation[-1]
            if "S" in firstEl:
                #strip off S and convert to int
                softstart = int(firstEl.strip("S"))
                if softstart > 38:
                    #trim sequence and quality score lines
                    sequence = sequence[:softstart]
                    quality = quality[:softstart]
                    #print to output
                    O.write("@{0}__{1}\n{2}\n+\n{3}\n".format(line[2], line[0], sequence, quality))
            if "S" in lastEl:
                softstop = int(lastEl.strip("S"))
                if softstop > 38:
                    sequence = sequence[-softstop:]
                    quality = quality[-softstop:]
                    #print to output
                    O.write("@{0}__{1}\n{2}\n+\n{3}\n".format(line[2], line[0], sequence, quality))
            fractionation = []
O.close()                    
                    
                
