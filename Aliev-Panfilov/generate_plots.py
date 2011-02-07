#!/usr/bin/python

import sqlite3, sys, re, os, math

def makePlotFile(f, testtype, datfiles):
    f.write("set terminal png size 1000, 700\n")
    f.write("set autoscale\n")
    f.write("unset log\n")
    f.write("unset label\n")
    f.write("set logscale x 2\n")
    f.write("set xtic auto\n")
    f.write("set ytic nomirror\n")
    f.write("set title \"Strong scaling runs\"\n")
    f.write("set xlabel \"Number of cores\"\n")
    f.write("set ylabel \"Speedup\"\n")
    f.write("plot ")
    for i in xrange(len(datfiles)):
        n,iters,df,srt = datfiles[i]
        if i > 0: f.write(",\\\n")
        f.write("\"%s\" using 1:(%f/$2) title \"n=%d, i=%d\" with linespoints" % (df,srt,n,iters))
    f.write("\n")    

def main(args):
    dbfile = "results.db"
    searchpath = "./"
    if len(args) > 0: dbfile = args[-1]

    # Open database
    db = sqlite3.connect(dbfile)
    c = db.cursor()

    # Generate plots from data. First get the list of test sets.
    tests = []
    c.execute("select distinct name, version from testresults")
    for n in c:
        tests.append( (n[0], n[1]) )

    # First generate plots for the strong scaling results. Put all strong scaling results in one plot. 
    for (test, version) in tests:
        c.execute("select problemsize, iters, threads, min(runningtime) from testresults where name='%s' and section='strong_scaling' and version='%s' group by problemsize, iters, threads order by problemsize" % (test, version))
        fnamePrefix = "strongscaling_" + test + "_" + version

        # Tuple of file name, file handle and serial execution time for each combination of n and i.
        files = {}
        for row in c:
            n = int(row[0])
            i = int(row[1])
            threads = int(row[2])
            rt = float(row[3])
            key = (n,i)
            if not key in files:
                fname = fnamePrefix + "_" + str(n) + "_" + str(i) + ".dat"
                files[key] = (fname, open(fname, "w"), 10e9)

            if threads == 1:
                files[key] = (files[key][0], files[key][1], rt)

            files[key][1].write(str(threads) + " " + str(rt) + "\n")

        plotfname = fnamePrefix + ".p"
        plotfile = open(plotfname, "w")

        datfiles = []
        for (n, i), (fname, f, srt) in files.iteritems():
            f.close()
            datfiles.append( (n,i,fname, srt) )
        makePlotFile(plotfile, "strong_scaling", datfiles)

            
    c.close()
    db.close()

if __name__ == "__main__":
    main(sys.argv[1:])
