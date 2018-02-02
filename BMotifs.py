import numpy as np
import pandas as pd

"""
Usage:
    bmotif counts occurrences of motifs in bipartite networks, as well as
    the number of times each node appears in each unique position within motifs.
    bmotif was intended for use in ecology but its methods are general and can be
    applied to any bipartite graph.

    To use, first:    
        >>> from BMotifs import Motifs

    Then, to generate all the matrices for the motif counts:
        >>> cm = Motifs(M=mat)
    where mat is a two-dimensional binary NumPy array. The nodes of the two
    bipartite layers are ordered along the rows and columns, respectively.
    
    To calculate the frequency with which each node in the matrix occurs in 
    different positions within motifs, use::
        >>> p=cm.pos_motifs(num_m)
        >>> p.position_motifs
    where num_m is a list containing the motif positions to count, between 1 and 148.

    p.position_motifs returns the resulting posistion motifs as a pandas table
    with a number of rows equal to the number of motif positions requested in
    num_m. The number of columns in the table is equal to the sum of the dimensions of
    mat (number of rows + number of columns). The first column (MotifID) gives
    the ID of the motif position, the Ax columns give the frequency with which 
    each row node occurs in each motif position, the Bx columns give the
    frequency with which each column node occurs in each motif position.
    NaN is returned when a node cannot occur in a given position because that 
    position can only be occupied by nodes in the other level. For example, 
    position 1 (MotifID 1) can only be occupied by column nodes and therefore 
    a count of NaN is returned for all row nodes for this position.
    
    To count the frequency with which motifs occurs in the matrix, use::
    >>> p=cm.motifs(num_m)
    >>> p.motifs
    where num_m is a list containing the motifs to count between 1 and 44.
    p.motifs return the resulting motifs as a pandas table with one row for
    each motif. The first column (MotifID) gives the ID of the motif,
    the second column gives the frequency with which that motif occurs in the network.
    
    The code is released under the MIT license.
        
"""

class Motifs(object):


    def __init__(self, M):
        """Initialize the parameters of the Motifs count.
        :Generating all the matrices for the motif counts
        :``M`` a two-dimensional binary NumPy array
          """
        print 'Iniziating the matrices ...'
        self.M = np.array(M, dtype=np.int64)
        self.check_input_selfrix_is_binary()
        [self.num_rows, self.num_columns] = self.M.shape
        self.motifs = None             #motifs
        self.position_motifs = None      #position motifs
        
        self.z=self.M.shape[0]
        self.znul=nans([self.z,1])
        self.p=self.M.shape[1]
        self.pnul=nans([self.p,1])
     
        self.jP=np.ones([self.p,1])
        self.jZ=np.ones([self.z,1])
        self.J=np.ones([self.z,self.p])
        self.JP=np.ones([self.p,self.p])
        self.JZ=np.ones([self.z,self.z])

        self.MT=self.M.transpose()
        self.N=self.J-self.M
        self.NT=self.N.transpose()


        self.dZ=np.matmul(self.M,self.jP)
        self.dP=np.matmul(self.MT,self.jZ)

        self.Z=np.matmul(self.M,self.MT)
        self.Y=np.matmul(self.M,self.NT)
        self.X=np.matmul(self.N,self.MT)

        self.P=np.matmul(self.MT,self.M)
        self.Q=np.matmul(self.MT,self.N)
        self.R=np.matmul(self.NT,self.M)
        
        self.JP3=np.ones([self.p,self.p,self.p]);
        self.JZ3=np.ones([self.z,self.z,self.z]);
        
        self.KP3=self.JP3;
        for i in range(0,self.p):
            for j in range(0,self.p):
                self.KP3[i,j,j]=0
                self.KP3[j,i,j]=0
                self.KP3[j,j,i]=0

        self.KZ3=self.JZ3
        for i in range(0,self.z):
            for j in range(0,self.z):
                self.KZ3[i,j,j]=0
                self.KZ3[j,i,j]=0
                self.KZ3[j,j,i]=0
                

        self.AZ=tensor_make(self.M,self.M);
        self.BZ=tensor_make(self.M,self.N);
        self.CZ=tensor_make(self.N,self.M);
        self.DZ=tensor_make(self.N,self.N);

        self.AP=tensor_make(self.MT,self.MT);
        self.BP=tensor_make(self.MT,self.NT);
        self.CP=tensor_make(self.NT,self.MT);
        self.DP=tensor_make(self.NT,self.NT);
        
        self.MTA = tensorR(self.MT,self.AZ);
        self.MTB = tensorR(self.MT,self.BZ);
        self.MTC = tensorR(self.MT,self.CZ);
        self.MTD = tensorR(self.MT,self.DZ);

        print 'Almonst Done ...'
        self.MA = tensorR(self.M,self.AP);
        self.MB = tensorR(self.M,self.BP);
        self.MC = tensorR(self.M,self.CP);
        self.MD = tensorR(self.M,self.DP);
        
        self.NTA = tensorR(self.NT,self.AZ);
        self.NTB = tensorR(self.NT,self.BZ);
        self.NTC = tensorR(self.NT,self.CZ);

        self.Na = tensorR(self.N,self.AP);
        self.NB = tensorR(self.N,self.BP);
        self.NC = tensorR(self.N,self.CP);
        
        self.MTAK = self.MTA*self.KP3;
        self.MTBK = self.MTB*self.KP3;
        self.MTCK = self.MTC*self.KP3;
        self.MTDK = self.MTD*self.KP3;

        self.MAK = self.MA*self.KZ3;
        self.MBK = self.MB*self.KZ3;
        self.MCK = self.MC*self.KZ3;
        self.MDK = self.MD*self.KZ3;

        self.NTAK = self.NTA*self.KP3;
        self.NTBK = self.NTB*self.KP3;
        self.NTCK = self.NTC*self.KP3;

        self.NaK = self.Na*self.KZ3;
        self.NBK = self.NB*self.KZ3;
        self.NCK = self.NC*self.KZ3;
        print 'Done'

        
        
    def check_input_selfrix_is_binary(self):
        """Check that the input selfrix is binary, i.e. entries are 0 or 1.
        :raise AssertionError: raise an error if the input selfrix is not
            binary
        """
        assert np.all(np.logical_or(self.M == 0, self.M == 1)), \
            "Input selfrix is not binary."
            
            
    def count_pos_motifs(self,V_motifs):
        V_motifs = np.array(V_motifs, dtype=np.int64)
        if not [isinteger(k) for k in V_motifs]:
            print 'check the Position motifs requested'
            return
        bigger=[x for x in V_motifs if x >= 149]
        lower=[x for x in V_motifs if x < 1]
        if len(bigger)>0 or len(lower)>0:
            print 'check the Position motifs requested'
            return
        
        Names=[]
        Names.append('MotifID')

        for i in range(0,self.z):
            Names.append('A'+str(i+1))
        for i in range(0,self.p):
            Names.append('B'+str(i+1))
            
        Motifout=np.zeros([len(V_motifs),self.p+self.z+1])
      
        
        for k in range(0,len(V_motifs)):
            mot=self.pos_motifs(V_motifs[k])
            Motifout[k,:]=np.append(V_motifs[k],mot)
            
        Motifout=pd.DataFrame(Motifout, columns=Names)
        Motifout.MotifID= Motifout.MotifID.astype(int)
        self.position_motifs=Motifout
        return Motifout
    
    
                
    def count_motifs(self,T_motifs):
        T_motifs = np.array(T_motifs, dtype=np.int64)
        if not [isinteger(k) for k in T_motifs]:
            print 'check the Position motifs requested'
            return
        bigger=[x for x in T_motifs if x >= 45]
        lower=[x for x in T_motifs if x < 1]
        if len(bigger)>0 or len(lower)>0:
            print 'check the Motifs requested'
            return
        
        Names=[]
        Names.append('MotifID')

        Names.append('Motif')
            
        Motifout=np.zeros([len(T_motifs),2])
      
        
        for k in range(0,len(T_motifs)):
            mot=self.c_motifs(T_motifs[k])
            Motifout[k,:]=np.append(T_motifs[k],mot)
            
        Motifout=pd.DataFrame(Motifout, columns=Names)
        Motifout.MotifID= Motifout.MotifID.astype(int)
        Motifout.Motif= Motifout.Motif.astype(int)
        self.motifs=Motifout
        return Motifout
    
    
    
    
    
    def pos_motifs(self,num_m):
        Mot={
            1: (self.dP,'z'),
            2: (self.dZ,'p'),
            3: (np.matmul(self.P,self.jP) - self.dP,'z'),
            4: (np.multiply(self.dZ,(self.dZ - self.jZ))/2,'p'),
            5: (np.multiply(self.dP, (self.dP - self.jP)) / 2,'z'),
            6: (np.matmul(self.Z,self.jZ)- self.dZ,'p'),
            7: (np.multiply(np.multiply(self.dP ,(self.dP - self.jP)),(self.dP - 2*self.jP)) / 6,'z'),
            8: (np.matmul(self.M, np.multiply((self.dP - self.jP), (self.dP - 2 * self.jP))) / 2,'p'),
            9: (np.matmul(np.multiply(self.P,self.R),self.jP),'z'),
            10: (np.matmul(np.multiply(self.P,self.Q),self.jP),'z'),
            11: (np.matmul(np.multiply(self.X,self.Z),self.jZ),'p'),
            12: (np.matmul(np.multiply(self.Y,self.Z),self.jZ),'p'),
            13: (np.matmul(np.multiply(self.P,(self.P - self.JP)),self.jP) / 2 - np.multiply(self.dP,(self.dP - self.jP)) / 2,'z'),
            14: (np.matmul(np.multiply(self.Z,(self.Z - self.JZ)),self.jZ)/2 - np.multiply(self.dZ,(self.dZ - self.jZ))/ 2,'p'),
            15: (np.matmul(self.MT,np.multiply((self.dZ - self.jZ),(self.dZ - 2*self.jZ)))/ 2,'z'),
            16: (self.dZ*(self.dZ - self.jZ) * (self.dZ - 2 * self.jZ) / 6,'p'),
            17: (self.dP * (self.dP - self.jP) * (self.dP - 2 * self.jP) * (self.dP - 3 * self.jP) / 24,'z'),
            18: (np.matmul(self.M ,((self.dP - self.jP) * (self.dP - 2 * self.jP) * (self.dP - 3*self.jP))) / 6,'p'),
            19: (np.matmul((self.P * self.R * (self.R - self.JP)) , self.jP) / 2,'z'),
            20: (np.matmul((self.P * self.Q * (self.Q - self.JP)) , self.jP) / 2,'z'),
            21: (np.matmul(self.N * np.matmul(self.M , ((self.Q - self.JP) * self.P)) , self.jP),'p'),
            22: (np.matmul((self.M * np.matmul(self.M , (self.Q * (self.Q - self.JP)))) , self.jP) / 2,'p'),
            23: (np.matmul((self.P * self.Q * self.R),self.jP),'z'),
            24: (np.matmul((self.N * np.matmul(self.M,(self.P * self.R))) , self.jP),'p'),
            25: (np.matmul((self.M * np.matmul(self.M , (self.Q * self.R))) , self.jP) / 2,'p'),
            26: (np.matmul((self.P * (self.P - self.JP) * self.R) ,self.jP) / 2,'z'),
            27: (np.matmul((self.P * (self.P - self.JP) * self.Q) , self.jP) / 2,'z'),
            28: (np.matmul((self.N * np.matmul(self.M , (self.P * (self.P - self.JP)))),  self.jP) / 2,'p'),
            29: (np.matmul((self.M * np.matmul(self.M , (self.Q * (self.P - self.JP)))) , self.jP),'p'),
            30: (np.matmul((self.P * (self.P - self.JP) * (self.P - 2 * self.JP)) , self.jP) / 6 -  self.dP * (self.dP - self.jP) * (self.dP - 2 * self.jP) / 6,'z'),
            31: ((np.matmul(self.M * np.matmul(self.M , ((self.P - self.JP) * (self.P - 2 * self.JP))) , self.jP) - np.matmul(self.M , ((self.dP - self.jP) * (self.dP - 2 * self.jP)))) / 4,'p'),
            32: (np.matmul((self.NT * np.matmul(self.MT , ((self.Y - self.JZ) * self.Z))) , self.jZ),'z'),
            33: (np.matmul((self.MT * np.matmul(self.MT,(self.Y * (self.Y - self.JZ)))) , self.jZ) / 2,'z'),
            34: (np.matmul((self.Z * self.X * (self.X - self.JZ)) , self.jZ )/ 2,'p'),
            35: (np.matmul((self.Z * self.Y * (self.Y - self.JZ)) , self.jZ) / 2,'p'),
            36: (np.matmul((self.NT * np.matmul(self.MT , (self.Z*self.X))) , self.jZ),'z'),
            37: (np.matmul((self.MT * np.matmul(self.MT , (self.Y * self.X))) , self.jZ) / 2,'z'),
            38: (np.matmul((self.Z * self.Y * self.X) , self.jZ),'p'),
            39: (np.matmul((self.NT * np.matmul(self.MT , (self.Z * (self.Z - self.JZ)))) ,self.jZ) / 2,'z'),
            40: (np.matmul((self.MT * np.matmul(self.MT , (self.Y * (self.Z - self.JZ)))) , self.jZ),'z'),
            41: (np.matmul((self.Z * (self.Z - self.JZ) * self.X) , self.jZ )/ 2,'p'),
            42: (np.matmul((self.Z * (self.Z - self.JZ) * self.Y) , self.jZ )/ 2,'p'),
            43: (np.matmul((self.MT * np.matmul(self.MT , ((self.Z - self.JZ) * (self.Z - 2 * self.JZ)))) , self.jZ) / 4 -  np.matmul(self.MT , ((self.dZ - self.jZ) * (self.dZ - 2 * self.jZ))) / 4,'z'),
            44: ((np.matmul((self.Z * (self.Z - self.JZ) * (self.Z - 2 * self.JZ)) , self.jZ)  -  self.dZ * (self.dZ - self.jZ) * (self.dZ - 2 * self.jZ) )/ 6,'p'),
            45: (np.matmul(self.MT , ((self.dZ - self.jZ) * (self.dZ - 2 * self.jZ) * (self.dZ - 3 * self.jZ)) )/ 6,'z'),
            46: (self.dZ * (self.dZ - self.jZ) * (self.dZ - 2 * self.jZ) * (self.dZ - 3 * self.jZ) / 24,'p'),
            47: (self.dP*(self.dP - self.jP)*(self.dP - 2*self.jP)*(self.dP - 3*self.jP)*(self.dP - 4*self.jP) / 120,'z'),
            48: (np.matmul(self.M , ((self.dP - self.jP)*(self.dP - 2*self.jP)*(self.dP - 3*self.jP)*(self.dP - 4*self.jP))) / 24,'p'),
            49: (np.matmul((self.P*self.R*(self.R - self.JP)*(self.R - 2*self.JP)) , self.jP) / 6,'z'),
            50: (np.matmul((self.P*self.Q*(self.Q - self.JP)*(self.Q - 2*self.JP)) , self.jP) / 6,'z'),
            51: (np.matmul(self.N*np.matmul(self.M , (self.P*(self.Q - self.JP)*(self.Q - 2*self.JP))) ,self.jP )/ 2,'p'),
            52: (np.matmul((self.M*np.matmul(self.M , (self.Q*(self.Q - self.JP)*(self.Q - 2*self.JP)))) , self.jP) / 6,'p'),
            53: (np.matmul((self.P*self.Q*self.R*(self.R - self.JP)) , self.jP) / 2,'z'),
            54: (np.matmul((self.P*self.Q*(self.Q - self.JP)*self.R) , self.jP) / 2,'z'),
            55: (np.matmul((self.M*np.matmul(self.N , (self.P*self.Q*(self.Q - self.JP)))) , self.jP) / 2,'p'),
            56: (np.matmul((self.M*np.matmul(self.N , (self.P*self.Q*(self.R - self.JP)))) , self.jP),'p'),
            57: (np.matmul((self.M*np.matmul(self.M , (self.Q*(self.Q - self.JP)*self.R))) , self.jP) / 2,'p'),
            58: (np.matmul((self.P*(self.P - self.JP)*self.R*(self.R - self.JP)) ,self.jP) / 4,'z'),
            59: (np.matmul((self.P*(self.P - self.JP)*self.Q*(self.Q - self.JP)) , self.jP) / 4,'z'),
            60: (np.matmul((self.N*np.matmul(self.M , ((self.Q - self.JP)*self.P*(self.P - self.JP)))) , self.jP )/ 2,'p'),
            61: (np.matmul((self.M*np.matmul(self.M , ((self.P - self.JP)*self.Q*(self.Q - self.JP)))) , self.jP )/ 2,'p'),
            62: (np.matmul((self.P*(self.P - self.JP)*self.Q*self.R) , self.jP) / 2,'z'),
            63: (np.matmul((self.N*np.matmul(self.M ,(self.P*(self.P - self.JP)*self.R))) ,self.jP) / 2,'p'),
            64: (np.matmul((self.M*np.matmul(self.M , ((self.P - self.JP)*self.Q*self.R))) , self.jP) / 2,'p'),
            65: (np.matmul((self.P*(self.P - self.JP)*(self.P - 2*self.JP)*self.R) , self.jP) / 6,'z'),
            66: (np.matmul((self.P*(self.P - self.JP)*(self.P - 2*self.JP)*self.Q) , self.jP) / 6,'z'),
            67: (np.matmul((self.N*np.matmul(self.M , (self.P*(self.P - self.JP)*(self.P - 2*self.JP)))) , self.jP) / 6,'p'),
            68: (np.matmul((self.M*np.matmul(self.M , ((self.P - self.JP)*(self.P - 2*self.JP)*self.Q))) , self.jP) / 2,'p'),
            69: (np.matmul((self.P*(self.P - self.JP)*(self.P - 2*self.JP)*(self.P - 3*self.JP)) , self.jP) / 24 - self.dP*(self.dP - self.jP)*(self.dP - 2*self.jP)*(self.dP - 3*self.jP) / 24,'z'),
            70: (np.matmul((self.M*np.matmul(self.M , ((self.P - self.JP)*(self.P - 2*self.JP)*(self.P - 3*self.JP)))) , self.jP) / 12 - np.matmul(self.M , ((self.dP - self.jP)*(self.dP - 2*self.jP)*(self.dP - 3*self.jP))) / 12,'p'),
            71: (np.sum(np.sum((self.NTB * (self.NTB - self.JP3) * self.MTA * self.KP3),axis=1),axis=1)/2,'z'),
            72: (np.sum(np.sum((self.MTD * (self.MTD - self.JP3) * self.MTA * self.KP3),axis=1),axis=1) / 4,'z'),
            73: (np.sum(np.sum((self.NB * (self.NB - self.JZ3) * self.MA * self.KZ3),axis=1),axis=1) / 2,'p'),
            74: (np.sum(np.sum((self.MD * (self.MD - self.JZ3) * self.MA * self.KZ3),axis=1),axis=1) / 4,'p'),
            75: (np.sum(np.sum(self.MTD * self.MTA * self.NTB * self.KP3,axis=1),axis=1),'z'),
            76: (np.sum(np.sum(self.MTA * self.NTB * self.NTC * self.KP3,axis=1),axis=1) / 2,'z'),
            77: (np.sum(np.sum(self.MD * self.MB * self.MC * self.KZ3,axis=1),axis=1) / 2,'p'),
            78: (np.sum(np.sum(self.MB * self.NB * self.Na * self.KZ3,axis=1),axis=1),'p'),
            79: (np.sum(np.sum(self.MTB * self.MTC * self.MTD * self.KP3,axis=1),axis=1) / 2,'z'),
            80: (np.sum(np.sum(self.MTB * self.NTB * self.NTA * self.KP3,axis=1),axis=1),'z'),
            81: (np.sum(np.sum(self.MA * self.MD * self.NB * self.KZ3,axis=1),axis=1),'p'),
            82: (np.sum(np.sum(self.MA * self.NB * self.NC * self.KZ3,axis=1),axis=1) / 2,'p'),
            83: (np.sum(np.sum(self.MTB * self.NTA * self.NTC * self.KP3,axis=1),axis=1),'z'),
            84: (np.sum(np.sum(self.NTA * self.MTB * self.MTD * self.KP3,axis=1),axis=1),'z'),
            85: (np.sum(np.sum(self.MTB * self.MTC * self.NTC * self.KP3,axis=1),axis=1),'z'),
            86: (np.sum(np.sum(self.MB * self.Na * self.NC * self.KZ3,axis=1),axis=1),'p'),
            87: (np.sum(np.sum(self.MD * self.MB * self.Na * self.KZ3,axis=1),axis=1),'p'),
            88: (np.sum(np.sum(self.NB * self.MB * self.MC * self.KZ3,axis=1),axis=1),'p'),
            89: (np.sum(np.sum(self.MTA * self.NTA * self.NTB * self.KP3,axis=1),axis=1),'z'),
            90: (np.sum(np.sum(self.MTA * self.MTB * self.NTB * self.KP3,axis=1),axis=1),'z'),
            91: (np.sum(np.sum(self.MTA * self.MTB * self.MTD * self.KP3,axis=1),axis=1),'z'),
            92: (np.sum(np.sum(self.MA * self.Na * self.NB * self.KZ3,axis=1),axis=1),'p'),
            93: (np.sum(np.sum(self.MB * self.NB * self.MA * self.KZ3,axis=1),axis=1),'p'),
            94: (np.sum(np.sum(self.MA * self.MB * self.MD * self.KZ3,axis=1),axis=1),'p'),
            95: (np.sum(np.sum(self.MTB * self.NTA * (self.NTA - self.JP3) * self.KP3,axis=1),axis=1) / 2,'z'),
            96: (np.sum(np.sum(self.MTB * (self.MTB - self.JP3) * self.NTA * self.KP3,axis=1),axis=1) / 2,'z'),
            97: (np.sum(np.sum(self.MTB * self.MTC * (self.MTC - self.JP3) * self.KP3,axis=1),axis=1) / 2,'z'),
            98: (np.sum(np.sum(self.MD * self.MA * self.Na * self.KZ3,axis=1),axis=1) / 2,'p'),
            99: (np.sum(np.sum(self.MA * self.NB * self.MC * self.KZ3,axis=1),axis=1),'p'),
            100: (np.sum(np.sum(self.MTA * self.NTA * (self.NTA - self.JP3) * self.KP3,axis=1),axis=1) / 4,'z'),
            101: (np.sum(np.sum(self.MTA * self.MTB * (self.MTB - self.JP3) * self.KP3,axis=1),axis=1) / 2,'z'),
            102: (np.sum(np.sum(self.MA * (self.MA - self.JZ3) * self.NB * self.KZ3,axis=1),axis=1) / 2,'p'),
            103: (np.sum(np.sum(self.MD * self.MA * (self.MA - self.JZ3) * self.KZ3,axis=1),axis=1) / 4,'p'),
            104: (np.sum(np.sum(self.MTD * self.MTA * self.NTA * self.KP3,axis=1),axis=1) / 2,'z'),
            105: (np.sum(np.sum(self.MTB * self.MTA * self.NTC * self.KP3,axis=1),axis=1),'z'),
            106: (np.sum(np.sum(self.MB * self.Na * (self.Na - self.JZ3) * self.KZ3,axis=1),axis=1) / 2,'p'),
            107: (np.sum(np.sum(self.MB * (self.MB - self.JZ3) * self.Na * self.KZ3,axis=1),axis=1) / 2,'p'),
            108: (np.sum(np.sum(self.MB * (self.MB - self.JZ3) * self.MC * self.KZ3,axis=1),axis=1) / 2,'p'),
            109: (np.sum(np.sum(self.MTA * (self.MTA - self.JP3) * self.NTB * self.KP3,axis=1),axis=1) / 2,'z'),
            110: (np.sum(np.sum(self.MTA * (self.MTA - self.JP3) * self.MTD * self.KP3,axis=1),axis=1) / 4,'z'),
            111: (np.sum(np.sum(self.MA * self.Na * (self.Na - self.JZ3) * self.KZ3,axis=1),axis=1) / 4,'p'),
            112: (np.sum(np.sum(self.MA * self.MB * (self.MB - self.JZ3) * self.KZ3,axis=1),axis=1) / 2,'p'),
            113: (np.sum(np.sum(self.MTB * self.MTC * self.NTA * self.KP3,axis=1),axis=1) / 2,'z'),
            114: (np.sum(np.sum(self.MB * self.MC * self.Na * self.KZ3,axis=1),axis=1) / 2,'p'),
            115: (np.sum(np.sum(self.MTA * self.MTB * self.NTA * self.KP3,axis=1),axis=1),'z'),
            116: (np.sum(np.sum(self.MTA * self.MTB * self.MTC * self.KP3,axis=1),axis=1) / 2,'z'),
            117: (np.sum(np.sum(self.MA * self.Na * self.MB * self.KZ3,axis=1),axis=1),'p'),
            118: (np.sum(np.sum(self.MA * self.MB * self.MC * self.KZ3,axis=1),axis=1) / 2,'p'),
            119: (np.sum(np.sum(self.MTA * (self.MTA - self.JP3) * self.NTA * self.KP3,axis=1),axis=1) / 4,'z'),
            120: (np.sum(np.sum(self.MTA * (self.MTA - self.JP3) * self.MTB * self.KP3,axis=1),axis=1) / 2,'z'),
            121: (np.sum(np.sum(self.MA * (self.MA - self.JZ3) * self.Na * self.KZ3,axis=1),axis=1) / 4,'p'),
            122: (np.sum(np.sum(self.MB * self.MA * (self.MA - self.JZ3) * self.KZ3,axis=1),axis=1) / 2,'p'),
            123: (np.sum(np.sum(self.MTA * (self.MTA - self.JP3) * (self.MTA - 2 * self.JP3) * self.KP3,axis=1),axis=1) / 12,'z'),
            124: (np.sum(np.sum(self.MA * (self.MA - self.JZ3) * (self.MA - 2 * self.JZ3) * self.KZ3,axis=1),axis=1) / 12,'p'),
            125: (np.matmul((self.NT * np.matmul(self.MT , (self.Z * (self.Y - self.JZ) * (self.Y - 2 * self.JZ)))) , self.jZ) / 2,'z'),
            126: (np.matmul((self.MT * np.matmul(self.MT , (self.Y * (self.Y - self.JZ) * (self.Y - 2 * self.JZ)))),  self.jZ) / 6,'z'),
            127: (np.matmul((self.Z * self.X * (self.X - self.JZ) * (self.X - 2 * self.JZ)) , self.jZ) / 6,'p'),
            128: (np.matmul((self.Z * self.Y * (self.Y - self.JZ) * (self.Y - 2 * self.JZ)) , self.jZ) / 6,'p'),
            129: (np.matmul((self.MT * np.matmul(self.NT , (self.Z * self.Y * (self.Y - self.JZ)))) , self.jZ) / 2,'z'),
            130: (np.matmul((self.MT * np.matmul(self.NT , (self.Z * self.Y * (self.X - self.JZ)))) , self.jZ),'z'),
            131: (np.matmul((self.MT * np.matmul(self.MT , (self.Y * (self.Y - self.JZ) * self.X))) , self.jZ) / 2,'z'),
            132: (np.matmul((self.Z * self.Y * self.X * (self.X - self.JZ)) , self.jZ) / 2,'p'),
            133: (np.matmul((self.Z * self.Y * (self.Y - self.JZ) * self.X) , self.jZ) / 2,'p'),
            134: (np.matmul((self.NT * np.matmul(self.MT , ((self.Y - self.JZ) * self.Z * (self.Z - self.JZ)))) , self.jZ) / 2,'z'),
            135: (np.matmul((self.MT * np.matmul(self.MT , ((self.Z - self.JZ) * self.Y * (self.Y - self.JZ)))) , self.jZ) / 2,'z'),
            136: (np.matmul((self.Z * (self.Z - self.JZ) * self.X * (self.X - self.JZ)) , self.jZ) / 4,'p'),
            137: (np.matmul((self.Z * (self.Z - self.JZ) * self.Y * (self.Y - self.JZ)) , self.jZ) / 4,'p'),
            138: (np.matmul((self.NT * np.matmul(self.MT , (self.Z * (self.Z - self.JZ) * self.X))) , self.jZ) / 2,'z'),
            139: (np.matmul((self.MT * np.matmul(self.MT , ((self.Z - self.JZ) * self.Y * self.X))) , self.jZ) / 2,'z'),
            140: (np.matmul((self.Z * (self.Z - self.JZ) * self.Y * self.X) , self.jZ) / 2,'p'),
            141: (np.matmul((self.NT * np.matmul(self.MT , (self.Z * (self.Z - self.JZ) * (self.Z - 2 * self.JZ)))) , self.jZ) / 6,'z'),
            142: (np.matmul((self.MT * np.matmul(self.MT , ((self.Z - self.JZ) * (self.Z - 2 * self.JZ) * self.Y))) , self.jZ) / 2,'z'),
            143: (np.matmul((self.Z * (self.Z - self.JZ) * (self.Z - 2 * self.JZ) * self.X) , self.jZ) / 6,'p'),
            144: (np.matmul((self.Z * (self.Z - self.JZ) * (self.Z - 2 * self.JZ) * self.Y) , self.jZ) / 6,'p'),
            145: (np.matmul((self.MT * np.matmul(self.MT , ((self.Z - self.JZ) * (self.Z - 2 * self.JZ) * (self.Z - 3 * self.JZ)))) , self.jZ) / 12 - np.matmul(self.MT , ((self.dZ - self.jZ) * (self.dZ - 2 * self.jZ) * (self.dZ - 3 * self.jZ))) / 12,'z'),
            146: (np.matmul((self.Z * (self.Z - self.JZ) * (self.Z - 2 * self.JZ) * (self.Z - 3 * self.JZ)) , self.jZ) / 24 - self.dZ * (self.dZ - self.jZ) * (self.dZ - 2 * self.jZ) * (self.dZ - 3 * self.jZ) / 24,'p'),
            147: (np.matmul(self.MT , ((self.dZ - self.jZ) * (self.dZ - 2 * self.jZ) * (self.dZ - 3 * self.jZ) * (self.dZ - 4 * self.jZ))) / 24,'z'),
            148: (self.dZ * (self.dZ - self.jZ) * (self.dZ - 2 * self.jZ) * (self.dZ - 3 * self.jZ) * (self.dZ - 4 * self.jZ) / 120,'p'),
                      
        }
        if Mot[num_m][1]=='z':
            out=np.append(self.znul.transpose(),Mot[num_m][0].transpose())
        else:
            out=np.append(Mot[num_m][0].transpose(),self.pnul.transpose())
        return out
    
    def c_motifs(self,num_m):
        Mot={
            1: np.sum(np.sum(self.M)),
            2: np.sum(self.dZ * (self.dZ - self.jZ)) / 2,
            3: np.sum(self.dP * (self.dP - self.jP)) / 2,
            4: np.sum(self.dP * (self.dP - self.jP) * (self.dP - 2 * self.jP)) / 6,
            5: np.sum(np.sum(self.Z * self.Y)),
            6: np.sum(np.sum(self.Z * (self.Z - self.JZ))) / 4 - np.sum(np.sum(self.dZ * (self.dZ - self.jZ))) / 4,
            7: np.sum(np.sum(self.dZ * (self.dZ - self.jZ) * (self.dZ - 2 * self.jZ))) / 6,
            8: np.sum(np.sum(self.dP * (self.dP - self.jP) * (self.dP - 2 * self.jP) * (self.dP - 3 * self.jP))) / 24,
            9: np.sum(np.sum(self.P * self.Q * (self.Q - self.JP))) / 2,
            10: np.sum(np.sum(self.P * self.Q * self.R)) / 2,
            11: np.sum(np.sum(self.P * (self.P - self.JP) * self.Q)) / 2,
            12: np.sum(np.sum(self.P * (self.P - self.JP) * (self.P - 2 * self.JP))) / 12 - np.sum(np.sum(self.dP * (self.dP - self.jP) * (self.dP - 2 * self.jP))) / 12,
            13: np.sum(np.sum(self.Z * self.Y * (self.Y - self.JZ))) / 2,
            14: np.sum(np.sum(self.Z * self.Y * self.X)) / 2,
            15: np.sum(np.sum(self.Z * (self.Z - self.JZ) * self.Y)) / 2,
            16: np.sum(np.sum(self.Z * (self.Z - self.JZ) * (self.Z - 2 * self.JZ))) / 12 - np.sum(np.sum(self.dZ * (self.dZ - self.jZ) * (self.dZ - 2 * self.jZ))) / 12,
            17: np.sum(np.sum(self.dZ * (self.dZ - self.jZ) * (self.dZ - 2 * self.jZ) * (self.dZ - 3 * self.jZ))) / 24,
            18: np.sum(np.sum(self.dP * (self.dP - self.jP) * (self.dP - 2 * self.jP) * (self.dP - 3 * self.jP) * (self.dP - 4 * self.jP)) ) / 120,
            19: np.sum(np.sum(self.P * self.Q * (self.Q - self.JP) * (self.Q - 2 * self.JP)) ) / 6,
            20: np.sum(np.sum(self.Q * self.P * self.R * (self.Q - self.JP)) ) / 2,
            21: np.sum(np.sum(self.Q * (self.Q - self.JP) * self.P * (self.P - self.JP)) ) / 4,
            22: np.sum(np.sum(self.Q * self.P * (self.P - self.JP) * self.R) ) / 4,
            23: np.sum(np.sum(self.P * (self.P - self.JP) * (self.P - 2 * self.JP) * self.Q) ) / 6,
            24: np.sum(np.sum(self.P * (self.P - self.JP) * (self.P - 2 * self.JP) * (self.P - 3 * self.JP))) / 48 - np.sum(np.sum(self.dP * (self.dP - self.jP) * (self.dP - 2 * self.jP) * (self.dP - 3 * self.jP))) / 48,
            25: np.sum(np.sum(np.sum(self.MAK * self.NBK * (self.NBK - self.JZ3))) ) / 4,
            26:  np.sum(np.sum(np.sum(self.MDK * self.MBK * self.MCK))) / 2,
            27:  np.sum(np.sum(np.sum(self.NBK * self.NCK * self.MAK))) / 2,
            28:  np.sum(np.sum(np.sum(self.MBK * self.MCK * self.NCK))),
            29:  np.sum(np.sum(np.sum(self.MAK * self.MBK * self.MDK))),
            30:  np.sum(np.sum(np.sum(self.MAK * self.MBK * self.NCK)))/2,
            31:  np.sum(np.sum(np.sum(self.MAK * self.NBK *(self.MAK - self.JZ3))))/4,
            32:  np.sum(np.sum(np.sum(self.MBK *( self.MBK - self.JZ3)*self.MCK)))/2,
            33:  np.sum(np.sum(np.sum(self.MBK *( self.MBK - self.JZ3)*self.MAK)))/4,
            34:  np.sum(np.sum(np.sum(self.NaK * self.MBK * self.MCK)))/6,
            35:  np.sum(np.sum(np.sum(self.MAK * self.MBK * self.MCK)))/2,
            36:  np.sum(np.sum(np.sum(self.MAK *( self.MAK -self.JZ3) * self.MBK)))/4,
            37:  np.sum(np.sum(np.sum(self.MAK *( self.MAK -self.JZ3) * (self.MAK -2*self.JZ3))))/36,
            38:  np.sum(np.sum(np.sum(self.Z *self.Y*( self.Y -self.JZ) * (self.Y-2*self.JZ))))/6,
            39:  np.sum(np.sum(np.sum(self.Z *( self.Y -self.JZ) * self.X *self.Y)))/2,
            40:  np.sum(np.sum(np.sum(self.Z *( self.Z -self.JZ) * self.Y *(self.Y-self.JZ))))/4,
            41:  np.sum(np.sum(np.sum(self.Z *( self.Z -self.JZ) * self.X * self.Y)))/4,
            42:  np.sum(np.sum(np.sum(self.Z *( self.Z -self.JZ)* (self.Z - 2 * self.JZ) * self.Y)))/6,
            43: np.sum(np.sum(np.sum(self.Z * (self.Z - self.JZ) * (self.Z - 2 * self.JZ) * (self.Z - 3 * self.JZ)))) / 48 -   np.sum(np.sum(np.sum(self.dZ * (self.dZ - self.jZ) * (self.dZ - 2 * self.jZ) * (self.dZ - 3 * self.jZ)))) / 48,
            44:  np.sum(np.sum(np.sum(self.dZ *( self.dZ - self.jZ) * (self.dZ - 2 * self.jZ) * (self.dZ - 3 * self.jZ) * (self.dZ - 4 * self.jZ))))/120,
        }
        
        return Mot[num_m]
        
        
        
        
            
        
def tensor_make(Ma,Mb):

    x=Ma.shape[0];
    y=Ma.shape[1];
    z=Mb.shape[1];

    output=np.ones([x,y,z])
    for j in range(0,y):
        for k in range(0,z):
            output[:,j,k]=np.multiply(Ma[:,j],Mb[:,k]);
    return output
        
        
        
def nans(shape, dtype=float):
    a = np.empty(shape, dtype)
    a.fill(np.nan)
    return a      
        
def tensorR(A,B):


    dimA=A.shape;
    dimB=B.shape;

    alongA=1;
    alongB=0;

    seqA=range(0,dimA[0])
    AllA=False
    seqAx=seqA
    del seqAx[alongA]
    permA=seqAx + [alongA]
    dimAx=list(dimA)
    del dimAx[alongA]
    Ax=np.reshape(A,[np.prod(dimAx),np.prod(dimA[alongA])])

    
    seqB=range(0,dimB[0])
    AllB=False
    seqBx=seqB
    del seqBx[alongB]
    permB=seqBx + [alongB]
    dimBx=list(dimB)
    del dimBx[alongB]
    Bx=np.reshape(B,[np.prod(dimB[alongB]),np.prod(dimBx)])


    R=np.matmul(Ax,Bx)
    
    Rx=np.reshape(R,dimAx+dimBx)
    
    return Rx
   
def isinteger(x):
    return np.equal(np.mod(x, 1), 0)





