import os
import subprocess as sb
import sys
usage = '\n- CRISPResso2 Docker Container -\n\t- Possible commands: \n\t  CRISPResso\n\t  CRISPRessoBatch\n\t  CRISPRessoPooled\n\t  CRISPRessoWGS\n\t  CRISPRessoCompare\n\t  CRISPRessoPooledWGSCompare\n\t  License\n\t- this special early-release docker version should be run like this:\n\tdocker run -v ${PWD}:/DATA -w /DATA -i kclem/trojan2 CRISPResso -r1 fastq1.fq -a AAAATTT \n'

if len(sys.argv)==1:

    print usage
    sys.exit(1)

if sys.argv[1]=='CRISPResso':
    sb.call(["/opt/conda/bin/python", "/CRISPResso/CRISPResso.py"]+ sys.argv[2:])
elif sys.argv[1]=='CRISPRessoBatch':
    sb.call(["/opt/conda/bin/python", "/CRISPResso/CRISPRessoBatch.py"]+ sys.argv[2:])
elif sys.argv[1]=='CRISPRessoCompare':
    sb.call(["/opt/conda/bin/python", "/CRISPResso/CRISPRessoCompare.py"]+ sys.argv[2:])
elif sys.argv[1]=='CRISPRessoPooled':
    sb.call(["/opt/conda/bin/python", "/CRISPResso/CRISPRessoPooled.py"]+ sys.argv[2:])
elif sys.argv[1]=='CRISPRessoWGS':
    sb.call(["/opt/conda/bin/python", "/CRISPResso/CRISPRessoWGS.py"]+ sys.argv[2:])
elif sys.argv[1]=='CRISPRessoPooledWGSCompare':
    sb.call(["/opt/conda/bin/python", "/CRISPResso/CRISPRessoPooledWGSCompare.py"]+ sys.argv[2:])
elif sys.argv[1]=='License':
    with open("EULA_for_CRISPResso2.txt", 'r') as fin:
        print fin.read()
else:
    print usage
