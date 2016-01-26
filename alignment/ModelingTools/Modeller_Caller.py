import subprocess
class modeller_caller:
    def __init__(self):
        self.modeller_executable = "python"
        self.process = None
    def run(self,script):
        try:
            process = subprocess.Popen([self.modeller_executable,script])
        except:
            print "Modeller or Script not found!"
        return process.wait()
