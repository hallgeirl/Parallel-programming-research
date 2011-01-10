#!/usr/bin/python
import sys
import os
import re

def runTests(ns, ts, geoms, reps, results, precision):
    r = re.compile("Running Time: ([0-9]+\.[0-9]+)")
    
    for (x,y) in geoms:
        if precision == "float":
            precision = " float=1"
        else:
            precision = ""
        os.system("make clean > /dev/null")
        os.system("make bx=" + str(x) + " by=" + str(y) + " device=3" + precision + " > /dev/null")
        for n in ns:
            for t in ts:
                print "n = " + str(n) + ", t = " + str(t)
                print "Thread geometry: " + str(x) + "x" + str(y)

                for rep in xrange(reps):
                    print "  Test run #:\t" + str(rep+1)
                    p = os.popen("./apf -n " + str(n) + " -t " + str(t))
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
    testruns = 10
    
    geomset = [(64,2), (16, 4), (16,16), (64,4), (128,2), (256,2)]

    if precision == "float":
        filename += "_float_"
    elif precision == "double":
        filename += "_double_"
    
    # Basic test for different thread geometries
    print "Running thread geometry tests..."

    # First for n=1023 with the thread geometry of 16x16
    runTests([1023], [5], [(16, 16)], testruns, results, precision)
    plot_threadgeoms = open(filename + "q1.dat", "w")
    plot_threadgeoms.write("# Question 1. n=1023, t=5, 16x16\n")
    key = (1023, 5, 16, 16)
    plot_threadgeoms.write("0 " + str(results[key]) + "\n")
    plot_threadgeoms.close()

    # Then n=1024
    runTests([1024], [5], [(16, 16)], testruns, results, precision)
    plot_threadgeoms = open(filename + "q3.dat", "w")
    plot_threadgeoms.write("# Question 3. n=1024, t=5, 16x16\n")
    key = (1024, 5, 16, 16)
    plot_threadgeoms.write("0 " + str(results[key]) + "\n")
    plot_threadgeoms.close()

    # Test for different thread geometries. n = 1023
    print "Running tests with different thread geomset"
    runTests([1023], [5], geomset, testruns, results, precision)
    plot_inctime = open(filename + "q2_1023.dat", "w")
    plot_inctime.write("# Question 2. n = 1023\n")
    plot_inctime.write("# time running_time\n")
    for i in xrange(len(geomset)):
        key = (1023, 5, geomset[i][0], geomset[i][1])
        plot_inctime.write("\"" + str(geomset[i][0]) + "x" + str(geomset[i][1]) + "\""  + " " + str(results[key]) + "\n")
    plot_inctime.close()   

    # Test for different thread geometries. n = 1024
    print "Running tests with different thread geomset"
    runTests([1024], [5], geomset, testruns, results, precision)
    plot_inctime = open(filename + "q2_1024.dat", "w")
    plot_inctime.write("# Question 2. n = 1024\n")
    plot_inctime.write("# time running_time\n")
    for i in xrange(len(geomset)):
        key = (1024, 5, geomset[i][0], geomset[i][1])
        plot_inctime.write("\"" + str(geomset[i][0]) + "x" + str(geomset[i][1]) + "\""  + " " + str(results[key]) + "\n")
    plot_inctime.close()   

    print "Results: " + str(results)

    
if __name__ == "__main__":
    main(sys.argv[1:])
