import re
import os
class pdb:
#  string file, chains;
    def __init__(self, arq):
        self.nomedoarq = arq
        self.file = file(arq, 'r')
        self.FileInList = self.file.readlines()        
    def Nome(self):
        return(self.file.name)
    def NomedaEstrutura(self):
        PrimeiraLinha = self.FileInList[0]
        if "HEADER" in PrimeiraLinha:
            return PrimeiraLinha[62:66].lower()
        else:
            print("Excecao: Arquivo PDB fora do Padrao")
    def nomedoarquivo(self):
        return(self.nomedoarq)
    def chains(self):
        List = []
        for i in self.LinhasOndeTemEscritoChain():
            texto = self.FileInList[i]
            if "CHAIN:" in texto:
                ChainPositionInTheString = texto.find(":") + 2
                if ";" in texto:
                    EndPositionInTheString = texto.find(";")
                else:
                    EndPositionInTheString = len(texto)-1                
                for j in range(ChainPositionInTheString,  EndPositionInTheString):
                    if texto[j] not in [' ',','] and texto[j] not in List:
                        List.append(texto[j])            
               
            else:
                #excecao
                print("error")
        return(List)
    def LinhasOndeTemEscritoChain(self):
        Linhas = []
        Contador = 0
        for i in self.FileInList: 

            if "CHAIN:" in i:
                Linhas.append(Contador)
            Contador+=1 
        return(Linhas)
    
    def hetatm(self):
        for i in self.FileInList:
            if  "HETATM" in i:
                return(True)
        return(False)
        
    def hoh(self):
       if self.hetatm():
            for i in self.FileInList:
                if "HOH" in i:
                    return(True) 
            return(False)
       else: 
           return(False)
    
    def HetatomsInPDB(self):
        if self.hetatm():
            List = []
            for i in self.FileInList:
                if ("HETATM" in i) and not ("REMARK" in i):
                    #Line = filter(None, re.split(" ", i))
                    if not i[17:20].replace(" ", "") in List:
                        List.append(i[17:20].replace(" ", ""))
            return(List)
        else:
            print("excecao")
            
    def ChangeHetatoms(self, Cadeia, Hetoriginal, Hetnovo):
        pdbsaida = ''
        if self.hetatm():
            for linha in self.FileInList:
                if ("HETATM" in linha) and not ("REMARK" in linha) :
                    #Line = filter(None, re.split(" ", linha))
                    if  (linha[17:20].replace(" ", "") == Hetoriginal) and ( linha[21] in Cadeia):
                        if len(Hetnovo) == 2:
                            Hetnovo = Hetnovo + " "
                        pdbsaida = pdbsaida + linha.replace(Hetoriginal,  Hetnovo)
                    else:
                        pdbsaida =  pdbsaida +linha
                else:
                        pdbsaida =  pdbsaida +linha
            return(pdbsaida)
                            

    
    def const(self,  hetatms,  chains):
        #if chains in self.chains() and hetatms in self.HetatomsInPDB():
            pdbsaida = ''
            for i in self.FileInList:
                Line = filter(None, re.split(" ", i))

                if(Line[0] != "CONECT"):
                    if ((Line[0] == "COMPND") and (Line[2] == "CHAIN:")):
                        texto = "COMPND   3 CHAIN: "
                        for letra in range(0, len(chains)):
                            texto = texto + chains[letra]
                            if (letra != len(chains)-1):
                                texto = texto + ","                        
                        pdbsaida = pdbsaida + texto + ";" + (62-2*len(chains))*(" ") + "\n"
                        #pdbsaida = pdbsaida 

                    else:
                        if ((i.startswith("ATOM")) or (i.startswith("HETATM"))):
                                 if((i[21] in chains) and (i[17:20].replace(" ", "") in hetatms[i[21]])):
                                    # pdbsaida = pdbsaida + "\n" + i
                                    pdbsaida = pdbsaida  + i
                                 else:
                                     pass
                                 if((i.startswith("ATOM")) and (i[21] in chains)):
                                    pdbsaida = pdbsaida + i
                                 else:
                                    pass
                        else:
                                # pdbsaida = pdbsaida + "\n" + i
                                pdbsaida = pdbsaida + i
            return(pdbsaida)

#    def const(self,  hetatms,  chains):
#        #if chains in self.chains() and hetatms in self.HetatomsInPDB():
#            pdbsaida = ''
#            for i in self.FileInList:
#                Line = filter(None, re.split(" ", i))
##                    if(Line[1] == "COMPND") and (Line[3] == "CHAIN:"):
#                        #Mudar cadeias
#                if(Line[0] != "CONECT"):
#                    if ((Line[0] == "COMPND") and (Line[2] == "CHAIN:")):
#                        texto = "COMPND   3 CHAIN: "
#                        for letra in range(0, len(chains)):
#                            texto = texto + chains[letra]
#                            if (letra != len(chains)-1):
#                                texto = texto + ","                        
#                        pdbsaida = pdbsaida + texto + ";" + (62-2*len(chains))*(" ") + "\n"
#                        #pdbsaida = pdbsaida 
#
#                    else:
#                        if ((Line[0] == "ATOM") or (Line[0] ==  "HETATM")):
#                                 if((Line[4][0] in chains) and (Line[3] in hetatms[Line[4][0]])):
#                                    # pdbsaida = pdbsaida + "\n" + i
#                                    pdbsaida = pdbsaida  + i
#                                 else:
#                                     pass
#                                 if((Line[0] == "ATOM") and (Line[4][0] in chains)):
#                                    pdbsaida = pdbsaida + i
#                                 else:
#                                    pass
#                        else:
#                                # pdbsaida = pdbsaida + "\n" + i
#                                pdbsaida = pdbsaida + i
#            return(pdbsaida)

    def HetatomsInChain(self,  chain):
        if self.hetatm():
            List = []
            for i in self.FileInList:
                if ("HETATM" in i) and not ("REMARK" in i):
                    #Line = filter(None, re.split(" ", i))
                    if ((not i[17:20].replace(" ", "") in List) and (chain == i[21])):
                        List.append(i[17:20].replace(" ", ""))
            return(List)
        else:
            print("excecao")

   
#al = pdb('/joao/1uij.pdb')
#al = pdb("/home/joao/IC/PDB/1uij.pdb")
#al = pdb("/home/joao/IC/1IKN.pdb")
#print al.nomedoarquivo()
#print al.Nome()
#print al.LinhasOndeTemEscritoChain()
#print al.chains()
#al = pdb("/joao/Alpha/src/teste1/models/1b8p.pdb")
#print al.ChangeHetatoms("A",   "NAG", "CO")

#print "Cadeias: ",    al.chains()
#print al.NomedaEstrutura()
#print "Heteroatomos? " ,  al.hetatm()
#print  "Agua? " ,  al.hoh()
#print "Heteroatomos: " ,  al.HetatomsInPDB()
#print al.HetatomsInChain('D')
#x = file("/joao/teste.txt",  "w")
#x.write(al.ChangeHetatoms("F",   "HOH", "CO"))
#x.close()

