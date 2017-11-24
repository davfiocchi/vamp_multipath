# coding=utf-8

import sys
import os
import numpy as np
import beat_evaluation_toolbox as be


def getBeatsFromFile(filePath):
	beatFile = open(filePath,"r")
	beats = []

	for line in beatFile:
		line = line.replace(lineTermination, "")
		line = line.split("\t")[0]
		beatTime = line.split(" ")[0]
		if beatTime is "":
			beatTime = line.split(" ")[1]
		beats.append(float('%.3f'%(float(beatTime))))

	beatFile.close()
	return beats

'''
Read the beats file
'''
ourPath = os.path.join(os.getcwd(), sys.argv[1]+os.sep)
lineTermination = "\n"

for (dirpath, dirnames, filenames) in os.walk(ourPath):
	ourBeatsList = []
	
	for file in filenames:
		print 'Analyzing (our) file: ' + str(file)
		ourBeatsList.append(np.array(getBeatsFromFile(os.path.join(ourPath,file))))
		#ourBeatsList.append(getBeatsFromFile(os.path.join(ourPath,file)))
	break

print 'Numero di nostri file letti: ' + str(len(ourBeatsList))
	
'''
Read the db file
'''
for (dirpath, dirnames, filenames) in os.walk(os.getcwd()):
	theirBeatsList = []
	
	for file in filenames:
		if str(file).endswith('.txt'):
			print 'Analyzing (their) file: ' + str(file)
			theirBeatsList.append(np.array(getBeatsFromFile(os.path.join(os.getcwd(),file))))
			#theirBeatsList.append(getBeatsFromFile(os.path.join(os.getcwd(),file)))
	break

print 'Numero di loro file letti: ' + str(len(theirBeatsList))

'''
Db evaluation
Va fatto per ogni canzone perchè non sappiamo come dare un nome diverso al file dei beat
'''
R = be.evaluate_db(theirBeatsList,ourBeatsList,'all',True)