#!/usr/bin/python

import sqlite3, sys, re, os

def parseFile(filename):
    lines = []
    with open(filename, "r") as f:
        lines = f.readlines()

    name = "unnamed_set"
    version = "unknown_version"
    section = "no_section"
    n = -1
    i = -1
    # Results stored as a mapping from (section, n, i) to (list of times)
    results = {}

    for l in lines:
        elements = l.strip().split(" ")
        if elements[0] == "testset_name":
            name = elements[1]
        elif elements[0] == "testset_version":
            version = elements[1]
        elif elements[0] == "section":
            section = "".join(elements[1:])
        elif elements[0][0] == "n" or elements[0][0] == "i":
            for e in elements:
                params = e.split("=")
                if params[0] == "n": n = int(params[1])
                elif params[0] == "i": i = int(params[1])
        elif elements[0] != "Nodes:":
            testresult = float(l)
            key = (section, n, i)
            if not key in results:
                results[key] = []
            results[key].append(testresult)
    return results
    

def main(args):
    dbfile = "results.db"
    searchpath = "./"
    if len(args) > 0: searchpath = args[-1] 
    if len(args) > 1: dbfile = args[-2]

    # Mapping from threads -> (list of results)
    results = {}

    files = os.listdir(searchpath)
    pattern = re.compile("[a-zA-Z_]+([0-9]+).sh.o[0-9]+")

    for f in files:
        match = pattern.match(f)
        if match: 
            threads = match.groups(1)[0]
            if not threads in results:
                results[threads] = []
            results[threads].append(parseFile(os.path.join(searchpath, f)))

    # Open database
    db = sqlite3.connect(dbfile)
    cursor = db.cursor()
    
    for k,v in results.iteritems():
        print str(k) + str(": ") + str(v)
    
    db.close()

if __name__ == "__main__":
    main(sys.argv[1:])
