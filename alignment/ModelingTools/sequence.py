class Sequence:
    def __init__(self,sequence_file):
        self.sequence_file = self.load_file_into_memory(sequence_file)

    def structures_names(self):
            file_in_tuple = self.sequence_file.splitlines()
            structures_list = []
            for eachline in file_in_tuple:
                if eachline.startswith('>P1;'):
                    struct = eachline[4:len(eachline)]
                    structures_list.append(struct)
            return structures_list
    def getStrname(self,number):
        if len(self.structures_names()[number]) == 5:
            return self.structures_names()[number][0:-1]
        else:
            return self.structures_names()[number]
                    
    def load_file_into_memory(self,arq):
        try:
            arquivo = open(arq,"rb")
        except IOError:
            print "Erro ao abrir o arquivo"
        arquivo.seek(0,2)
        tamanho = arquivo.tell()
        arquivo.seek(0)
        buffer = arquivo.read(tamanho)
        arquivo.close()
        return buffer

