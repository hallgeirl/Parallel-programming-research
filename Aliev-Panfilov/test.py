#!/usr/bin/python
import sys
import os
import re

def runTests(ns, ts, geoms, reps, results, precision):
    r = re.compile("Running Time: ([0-9]+\.[0-9]+)")

    if precision == "float":
        precision = " float=1"
    else:
        precision = ""

    os.system("make clean > /dev/null")
    os.system("make > /dev/null 2> /dev/null")
    
    for (x,y) in geoms:
        for n in ns:
            for t in ts:
                print "n = " + str(n) + ", t = " + str(t)
                print "Block geometry: " + str(x) + "x" + str(y)

                for rep in xrange(reps):
                    print "  Test run #:\t" + str(rep+1)
                    p = os.popen("./apf -n " + str(n) + " -t " + str(t) + " -i " + str(x) + " -j " + str(y))
                    key = (n, t, x, y)
                    for line in p.readlines():
                        if r.match(line):
                            time = float(r.findall(line)[0])
                            print "  Time: \t" + str(time) + " seconds."
                            if (key in results and time < results[key]) or key not in results:
                                print "  <<New best time>>"
                                results[key] = time
                            print ""
        

def main(args):
    # Keys: (n, t, x, y)-tuple
    # Values: (runtime, normalized runtime)
    results = {}
    filename = args[0]
    precision = args[1]
    testruns = 3
    
    nts = [(100, 5500), (200, 1500), (500, 90), (1000, 6), (2000, 0.4)];
    if precision == "float":
        filename += "_float_"
    elif precision == "double":
        filename += "_double_"
    
    # Basic test for different thread geometries
    print "Running tests..."
    plot_inctime = open(filename + "plotdata.dat", "w")
    plot_inctime.write("# time running_time\n")
    for threads in [1,2,4]:
        os.putenv("OMP_NUM_THREADS", str(threads))
        print str(threads) + " threads"
        plot_inctime.write("# " + str(threads) + " threads")
        for (n,t) in nts:
            x = n+1
            y = 30
            runTests([n], [t], [(x, y)], testruns, results, precision)
            plot_inctime.write("# n=" + str(n) + " t=" + str(t) + "\n")
            key = (n, t, x, y)
            plot_inctime.write("\"" + str(x) + "x" + str(y) + "\""  + " " + str(results[key]) + "\n")
    plot_inctime.close()   

    print "Results: " + str(results)

    
if __name__ == "__main__":
    main(sys.argv[1:])
