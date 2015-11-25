#!/usr/bin/python

import os, re, shutil

def run (cmd):
	print (cmd)
	os.system (cmd)

datadir = "RHF_example_data"
outputdir = "test"

shutil.rmtree(outputdir)
os.mkdir (outputdir)

result = open(outputdir+'/result.html', 'w')
result.write ("<html><body>\n")

for k in ["0","2"]:	
	for subdir, dirs, files in os.walk(datadir):
		for file in files:
			if re.search ('_hist.exr', file) :
				basename = re.match ('(.*)_hist.exr', file).group(1)
				input = datadir + "/" + basename + ".exr"
				output = outputdir + "/" + basename + "_k"+k + "_denoised.exr"
				outputpng = outputdir + "/" + basename + "_k"+k + "_denoised.png"
				diff = outputdir + "/" + basename + "_k"+k + "_diff.exr"
				diffpng = outputdir + "/" + basename + "_k"+k + "_diff.png"
				ref = datadir + "/" + basename + "_k"+k + "_denoised.exr"

				# Denoise
				run ("./rhf -k " + k + " -d 0.8 -h " + datadir+"/"+file + " " + input + " " + output)
				run ("./exrtopng " + output + " " + outputpng)

				# Difference
				run ("./exrdiff " + output + " " + ref + " 0.1 " + diff)
				run ("./exrtopng " + diff + " " + diffpng)

				result.write ("<h1>" + basename + " k=" + k + "</h1>\n")
				result.write ("<p><img src='../"+outputpng+"'> <img src='../"+diffpng+"'></p>\n")

result.write ("</body><html>\n")
result.close ()

os.system ("gnome-open test/result.html&")