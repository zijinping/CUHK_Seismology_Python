import struct
import numpy as np

class sac_file:
    def __init__(self,sFile):
        f = open(sFile,'rb')
        hdrBin = f.read(632)
        sfmt='f'*70+'I '*40+'8s '*22+'16s'
        hdrFmt=struct.Struct(sfmt)
        self.m_header = hdrFmt.unpack(hdrBin)
        self.delta=round(self.m_header[0],3)
        self.b=self.m_header[5]
        
        self.e=self.m_header[6]
        npts=int(self.m_header[79])
        

        fmt_data='f'*npts
        dataFmt=struct.Struct(fmt_data)
        dataBin=f.read(4*npts)
        f.close()
        self.data=dataFmt.unpack(dataBin)
        self.data=np.array(self.data)
        self.data.dtype='float32'

