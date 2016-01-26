import os
class RunCommand:
    def __init__(self, parameters):
        # self.command = [command]
        self.parameters = parameters
        self.process = None
        # fullcommand = self.modeller_executable  + self.parameters
    def run(self):
        # try:
            
            os.popen(" ".join(self.parameters))
            # self.process = subprocess.Popen(self.parameters)
        # except:
        #     print "Modeller or Script not found!"
        # return self.process.wait()

