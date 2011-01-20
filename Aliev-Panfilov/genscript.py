#!/usr/bin/python
#Generates run.sh scripts.
import sys
import os
import re

def main(args):
    testruns = 3
    
    # Create output directory if it doesn't exist
    if not os.path.exists("output"):
        os.mkdir("output")
        
    # Common header. All scripts use the same.
    header = "#!/bin/bash\n#PBS -A csd156\n#PBS -q large\n#PBS -l walltime=00:15:00\n"
    
    # Tails (that comes after the #PBS statements)
    common_tail = "cd $PBS_O_WORKDIR\nexport LD_PRELOAD=/opt/ScaleMP/libvsmpclib/0.1/lib64/libvsmpclib.so\nexport PATH=/opt/ScaleMP/numabind/bin:$PATH\n"
    omp_tail = "export KMP_AFFINITY=compact,verbose,0,`numabind --offset=16`\n"
    pthreads_tail = "cpu_start=`numabind --offset=$num_threads`\ncpu_end=`echo \"$cpu_start + $num_threads - 1\" | bc`\ntaskset -c $cpu_start-$cpu_end ./pthread.ex\n"

    #nts = [(100, 5500), (200, 1500), (500, 90), (1000, 6), (2000, 0.4), (5000, 0.01)];
    nts = [(1000, 6), (2000, 0.4), (5000, 0.01)];

    # Generate scripts for both pthreads and openmp
    for omp in [True, False]:
        for threads in [1,2,4,8,16,32]:
            filename = "run_"
            if omp:
                filename += "openmp"
            else:
                filename += "pthreads"
            filename += "_" + str(threads) + ".sh"
            f = open(filename, "w")
            f.write(header)
            f.write("#PBS -l nodes=1:ppn=" + str(threads) + "\n")
            if omp:
                f.write(omp_tail)
                f.write("export OMP_NUM_THREADS=" + str(threads) + "\n")
            else:
                f.write(pthreads_tail)
                f.write("num_threads=16\n")
                
            f.write(common_tail)
            
            # Jobs
            for n,t in nts:
                for foo in xrange(testruns):
                    f.write("./apf -n " + str(n) + " -t " + str(t))
                    if omp:
                        f.write(" -i " + str(n+1) + " -j 30 ")
                    else:
                        f.write(" -x 1 -y " + str(threads))
                    f.write("|grep Running|sed \"s/Running Time: //g\"|sed \"s/ sec.//g\"\n")
                        
if __name__ == "__main__":
    main(sys.argv[1:])
