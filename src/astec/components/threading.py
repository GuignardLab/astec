
import os
import subprocess
import threading
import time

class CompressFile(threading.Thread):
    def __init__(self, uncompressedfilename):
        threading.Thread.__init__(self)
        self._uncompressedfilename = uncompressedfilename

    def run(self):
        print("       .. compressing '" + os.path.basename(self._uncompressedfilename) + "'")
        command_line = "gzip " + str(self._uncompressedfilename)
        subprocess.call(command_line, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        print("       .. compressing of '" + os.path.basename(self._uncompressedfilename) + "' is done")


def waitForRunningThreadToStop(maxthread=5):
    # there is also the man thread which is running
    while threading.active_count() > maxthread:
        print("       .. there are " + str(threading.active_count()) + " active threads")
        time.sleep(10)
    return