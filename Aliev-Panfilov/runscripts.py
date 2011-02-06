#!/usr/bin/python

import os, subprocess, re, time

def machineIsFree():
	proc = subprocess.Popen("qstat|grep -v \"C large\"|grep large", shell=True, stdout=subprocess.PIPE)
	output = proc.stdout.read()
	proc.wait()
	if len(output.strip()) == 0:
		return True
	return False

def runJob(jobfile):
	proc = subprocess.Popen("qsub " + jobfile, shell=True, stdout=subprocess.PIPE)
	proc.wait()

def main():
	files = os.listdir("./")
	scripts = []
	for f in files:
		if f[-3:] == ".sh":
			scripts.append(f)
	os.system("clear")	
	print "Found scripts " + str(scripts)
	for s in scripts:
		print "Waiting for machine to be free..."
		while not machineIsFree():
			time.sleep(5)
		print "Running job " + s
		runJob(s)
		time.sleep(5)

		

if __name__=="__main__":
	main()
