from __future__ import print_function
import glob
import os

filelist = glob.glob("8*.hgh")

# LDA case is more tricky, some simply needs manual renaming

XC = "LDA"
for fil in filelist:
    nosymbol = fil.split(".")[0]
    print("")
    print("Filename = ", fil)
    if len(nosymbol) == 2:
        s = nosymbol[1]
        symbol = s[:1].upper() + s[1:]
    elif len(nosymbol) == 3:
        symbol = nosymbol[1]
        try:
			      s = int(symbol)
        except:
            s = 1
        if s > 0:
            symbol = nosymbol[2].upper()
        else:
            symbol == "X"
    elif len(nosymbol) == 4:
        s = nosymbol[2] + nosymbol[3]
        symbol = s[:1].upper() + s[1:]

    print("symbol = ", symbol)
    charge = fil.split(".")[1]
    newname = symbol + "." + charge + ".LDA.hgh"
    print("newname = ", newname)
    os.rename(fil,newname)

#XC = "GGA"
#for fil in filelist:
#    symbol = fil.split("-")[0]
#    charge = fil.split("-")[1].split(".")[0].replace("q","")
#    newname = symbol + "." + charge + ".GGA.hgh"
#    print("symbol = ", symbol, " charge = ", charge)
#    os.rename(fil,newname)


