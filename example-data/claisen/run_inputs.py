import os
from os import path

for FileName in sorted(list(os.listdir('.'))):
   BaseName, Ext = path.splitext(FileName)
   if Ext != '.inp':
      continue
   Cmd = "$MOLDEV2/bin/molpro -n 4 -W /tmp -d /tmp %s.inp" % BaseName
   print "!%s" % Cmd
   os.system(Cmd)
