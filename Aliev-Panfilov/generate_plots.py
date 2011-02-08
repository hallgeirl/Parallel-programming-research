#!/usr/bin/python

import sqlite3, sys, re, os, math

def makePlotFile(f, outputtype, datfiles, title, output_filename):
    f.write("set terminal png size 1000, 700\n")
    f.write("set output \"%s\"\n" % output_filename)
    f.write("set autoscale\n")
    f.write("unset log\n")
    f.write("unset label\n")
    f.write("set logscale x 2\n")
    f.write("set xtic auto\n")
    f.write("set ytic nomirror\n")
    f.write("set title \"" + title + "\"\n")
    f.write("set xlabel \"Number of cores\"\n")
    if outputtype == "speedup":
        f.write("set ylabel \"Speedup\"\n")
    elif outputtype == "efficiency":
        f.write("set ylabel \"Parallel efficiency\"\n")
    else:
        f.write("set ylabel \"Running time\"\n")
        f.write("set logscale y 2\n")
    
    f.write("plot ")
    for i in xrange(len(datfiles)):
        n,iters,df,srt = datfiles[i]
        if i > 0: f.write(",\\\n")
        f.write("\"%s\" using 1:" % (df))
        if outputtype == "speedup":
            f.write("(%f/$2)" % (srt))
        elif outputtype == "efficiency":
            f.write("((%f/$2)/$1)" % (srt));
        else:
            f.write("2")
        f.write(" title \"n=%d, i=%d\" with linespoints" % (n,iters))
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
        for (section, probsize) in [("strong_scaling", "problemsize"), ("weak_scaling", "problemsize_perthread")]:
            c.execute("select %s, iters, threads, min(runningtime) from testresults where name='%s' and section='%s' and version='%s' group by %s, iters, threads order by %s" % (probsize, test, section, version, probsize, probsize))
            fnamePrefix = test + "_" + section + "_" + version

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

            for outputtype in ["efficiency", "speedup", "runningtime"]:
                plotfname = fnamePrefix + "_" + outputtype + ".p"
                plotfile = open(plotfname, "w")

                datfiles = []
                for (n, i), (fname, f, srt) in files.iteritems():
                    f.close()
                    datfiles.append( (n,i,fname, srt) )
                title = ""
                datfiles = sorted(datfiles, key=lambda x: int(x[0]))
                if section == "strong_scaling": title = "Strong scaling runs."
                elif section == "weak_scaling": title = "Weak scaling runs. Problem size n is per thread."
                title += " Version: " + version
                makePlotFile(plotfile, outputtype, datfiles, title, plotfname[:-2] + ".png")
                plotfile.close()
                os.system("gnuplot " + plotfname) 
            
    c.close()
    db.close()

if __name__ == "__main__":
    main(sys.argv[1:])
