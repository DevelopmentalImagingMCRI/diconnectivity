import subprocess
import re

def getDIConnEnv():
# read environment variables
	p = subprocess.Popen(['which', 'DIConnEnv.sh'], stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	(stdoutdata, stderrdata) = p.communicate()

	if stdoutdata == None:
		return None
	else:
		stdoutdata = stdoutdata.strip()
		fp = open(stdoutdata, 'r')
		environmentVariables = dict()

		for curLine in fp.xreadlines():
			mat = re.match('^\s*(\w+)=(\w+)$', curLine)
			if mat != None:
				environmentVariables[mat.group(1)] = mat.group(2)
		
		fp.close()
		return environmentVariables

