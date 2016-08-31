def list_unique(inlist):
    outlist = []
    for x in inlist:
        if not x in outlist:
            outlist.append(x)
    return outlist
