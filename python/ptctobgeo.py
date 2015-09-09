import partio
import os
fileLocation = os.path.dirname(__file__)
print fileLocation
for i in range (0,100):
    p=partio.read(fileLocation+"SnowF"+str(i).zfill(4)+".ptc")
    partio.write(fileLocation+"SnowF"+str(i).zfill(4)+".bgeo",p)
