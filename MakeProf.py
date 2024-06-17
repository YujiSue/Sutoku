import os
import sys
import json
prvheader = ''
with open(sys.argv[1]) as f:
  prvheader = f.read()
appprof = json.load(open(sys.argv[2]))
header = prvheader.replace('%_PROFILE_%', json.dumps(appprof).replace('"','\\"'))
with open(sys.argv[3], "w") as f:
  f.write(header)
print("completed.")