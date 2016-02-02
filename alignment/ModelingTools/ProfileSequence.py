import re

class ProfileSequence():
    def __init__(self, estrutura):
        self.estrutura = estrutura
        self.estruturaparticionada = filter(None,  re.split(" ",  self.estrutura))
    def __str__(self):
        return float(self.estruturaparticionada[10])
    def name(self):
        return self.estruturaparticionada[1]
    def type(self):
        return self.estruturaparticionada[2]
    def startofStructure(self):
        return self.estruturaparticionada[5]
    def endofStructure(self):
        return self.estruturaparticionada[6]
    def sequenceSize(self):
        return self.estruturaparticionada[9]
    def identity(self):
        return self.estruturaparticionada[10]
    def sequence(self):
        return self.estruturaparticionada[12]