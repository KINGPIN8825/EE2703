'''EE2703 - assignment 1 - EE20B056 KATHIR PAGALAVAN - 26/01/2022'''

from sys import argv,exit

if len(argv)!=2:	#CHECKING IF THE NUMBER OF ARGUMENTS SUPPLIED IS CORRECT
	print('\nUsage :python3 %s <input-file>'% argv[0])
	exit()

CIRCUIT = '.circuit'	#DECLARING CONSTANTS TO COMPARE
END = '.end'

try:
	with open(argv[1]) as f:				#FILE OPENING
		lines = f.readlines()
		start=-1;end=0
		for line in lines:
			if CIRCUIT == line[:len(CIRCUIT)]:	#CHECKING FOR ".circuit"
				start=lines.index(line)	
			if END == line[:len(END)]:		#CHECKING FOR ".end"
				end=lines.index(line)
		if end>start and start>=0:			#CHECKING FOR VALID FORMAT
			for line in reversed(lines[start+1:end]):
				c=reversed(line.split('#')[0].split())	
				c = "  ".join(c)
				print(c)
				
			print("\nELEMENTAL ANALYSIS:\n")
			for i in range(start+1, end):
				line = lines[i]
				line=line.split("#")[0]
				tokens=line.split()
				if tokens[0][0]=='R':	#IF-TREE TO DETECT THE KIND OF ELEMENT
					print("Resistor of value "+tokens[3]+"ohms from "+tokens[1]+" to "+tokens[2]+"\n")
				elif tokens[0][0]=='L':
					print("Inductor of value "+tokens[3]+"henries from "+tokens[1]+" to "+tokens[2]+"\n")
				elif tokens[0][0]=='C':
					print("Capacitor of value "+tokens[3]+"farads from "+tokens[1]+" to "+tokens[2]+"\n")
				elif tokens[0][0]=='V':
					print("Voltage source of value "+tokens[3]+"volts from "+tokens[1]+" to "+tokens[2]+"\n")
				elif tokens[0][0]=='I':
					print("Current source of value "+tokens[3]+"amperes from "+tokens[1]+" to "+tokens[2]+"\n")
				elif tokens[0][0]=='E':
					print("Voltage source depending on voltage between "+tokens[3]+" and "+tokens[4]+"from "+tokens[1]+" to "+tokens[2]+"with value "+tokens[5]+"\n")
				elif tokens[0][0]=='G':
					print("Current source depending on voltage between "+tokens[3]+" and "+tokens[4]+"from "+tokens[1]+" to "+tokens[2]+"with value "+tokens[5]+"\n")
				elif tokens[0][0]=='H':
					print("Voltage source depending on current through the voltage source "+ tokens[3]+"with value "+tokens[4]+" from "+ tokens[1]+"to"+tokens[2]+"\n")
				elif tokens[0][0]=='I':
					print("Current source depending on current through the voltage source "+ tokens[3]+"with value "+tokens[4]+" from "+ tokens[1]+"to"+tokens[2]+"\n"	)
				else:
					print("Invalid element"+"\n")
		else:
			print("Invalid circuit definition")
			exit()

except FileNotFoundError:	#EXCEPTION TO HANDLE INEXISTENCE OF FILE WITH THE GIVEN NAME
	print("Invalid File.")
	exit()
			
