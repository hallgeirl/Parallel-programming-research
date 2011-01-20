#!/usr/bin/python
#Generates run.sh scripts.
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
                print "Thread geometry: " + str(x) + "x" + str(y)

                for rep in xrange(reps):
                    print "  Test run #:\t" + str(rep+1)
                    p = os.popen("./apf -n " + str(n) + " -t " + str(t) + " -x " + str(x) + " -y " + str(y))
                    key = (n, t, x, y)
                    found = False
                    lines = p.readlines()
                    for line in lines:
                        if r.match(line):
                            found = True
                            time = float(r.findall(line)[0])
                            print "  Time: \t" + str(time) + " seconds."
                            if (key in results and time < results[key]) or key not in results:
                                print "  <<New best time>>"
                                results[key] = time
                            print ""
                    if not found:
                        print "Some error occured. Program output: "
                        print lines
        

def main(args):
    testruns = 3
    
    # Create output directory if it doesn't exist
    if not os.path.exists("output"):
        os.mkdir("output")
        
    # Common header. All scripts use the same.
    header = "#!/bin/bash\n#PBS -A csd156\n#PBS -q large\n#PBS -o output/$PBS_JOBNAME_$PBS_JOBID\n"
    
    # Tails (that comes after the #PBS statements)
    common_tail = "cd $PBS_O_WORKDIR\nexport LD_PRELOAD=/opt/ScaleMP/libvsmpclib/0.1/lib64/libvsmpclib.so\nexport PATH=/opt/ScaleMP/numabind/bin:$PATH\n"
    omp_tail = "export KMP_AFFINITY=compact,verbose,0,`numabind --offset=16`\n"
    pthreads_tail = "cpu_start=`numabind --offset=$num_threads`\ncpu_end=`echo \"$cpu_start + $num_threads - 1\" | bc`\ntaskset -c $cpu_start-$cpu_end ./pthread.ex\n"

    #nts = [(100, 5500), (200, 1500), (500, 90), (1000, 6), (2000, 0.4), (5000, 0.01)];
    nts = [(1000, 6), (2000, 0.4), (5000, 0.01)];

    # Generate scripts for both pthreads and openmp
    for omp in [True, False]:
        for threads in [1,2,4,8,16,32]:
            f = open("run_" + ("openmp" if omp else "pthreads") + "_" + str(threads) + ".sh", "w")
            f.write(header)
            f.write("#PBS -l nodes=1:ppn=" + str(threads) + "\n")
            f.write(omp_tail if omp else pthreads_tail)
            f.write(("export OMP_NUM_THREADS=" + str(threads) + "\n") if omp else "num_threads=16\n")
            f.write(common_tail)
            
            # Jobs
            for n,t in nts:
                for foo in xrange(testruns):
                    f.write("./apf -n " + str(n) + " -t " + str(t))
                    if omp:
                        f.write(" -i " + str(n+1) + " -j 30 ")
                    else:
                        f.write(" -x 1 -y " + str(threads))
                    f.write("\n")
                        
if __name__ == "__main__":
    main(sys.argv[1:])
