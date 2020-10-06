import time
import sys

for i in range(100):
    time.sleep(0.1)
    sys.stdout.write("\r%02d%%" % (i+1))
    sys.stdout.flush()

