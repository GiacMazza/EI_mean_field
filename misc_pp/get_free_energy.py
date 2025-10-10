#!/usr/bin/python3
from  numpy import *
import numpy as np
import argparse
# from scipy.interpolate import interp1d
import sys
import scipy



#

#matplotlib.font_manager.findSystemFonts(fontpaths=None, fontext='ttf')
#[f.name for f in matplotlib.font_manager.fontManager.ttflist]

#

parser=argparse.ArgumentParser(description='plots variables')
parser.add_argument('--figNAME',dest='figNAME',default='fig',help='produce a pdf output')
parser.add_argument('--outPDF',dest='outPDF',action='store_true',help='produce a pdf output')
parser.add_argument('--outSVG',dest='outSVG',action='store_true',help='produce a SVG output')
args = parser.parse_args()




#
h_evps=4.135667662e-03 # Planck in eV*ps
Wev=0.5              # Bandwidth in eV
Wband=2.0


DIRPLOT='/home/giacomo/DATA/TNS_ei/TRS_2025/devel/long_range/driver_tests/PHN_tests/with_Bfield/UTa0.50_Uni1.50_V1.00_sweep_TEMP/gphn_0.000/read_cdsb/integrate_energy/'
PLOTFILE=DIRPLOT+'ene_VS_temp'
print (PLOTFILE)
TK,EeV = loadtxt(PLOTFILE,usecols=(0,1),unpack=True)

#interpolate on a log scale
log_TK=np.logspace(0,np.log10(TK[size(TK)-1]),10000)
# OUTFILE=DIRPLOT+'log_list_temp'
# np.savetxt(OUTFILE, np.c_[log_ene])

#xnew = np.linspace(0, 10, num=1001) #
splE = np.interp(log_TK,TK, EeV)



kb=8.61733326E-5
omega0=0.01
nph=1./(exp(omega0/kb/TK)-1.0)
#EeV=omega0*nph
#
#
#

ET2=EeV/TK/TK
ene_int=0.0*ET2
cum_int=0.0
for i in range(size(ET2)-1):
    cum_int=cum_int+0.5*(ET2[i+1]+ET2[i])*(TK[i+1]-TK[i])
    ene_int[i+1]=cum_int
    print (i)


OUTFILE=DIRPLOT+'int_ene_over_t2'
np.savetxt(OUTFILE, np.c_[TK,ene_int,ET2])
free_energy=TK/TK[0]*EeV[0]-TK*ene_int
entropy=(EeV-free_energy)/TK
check_ent=kb*((1+nph)*log(1+nph)-nph*log(nph))

OUTFILE=DIRPLOT+'free_ene_VS_temp'
np.savetxt(OUTFILE, np.c_[TK,free_energy,entropy,check_ent])



splET2=splE/log_TK/log_TK
spl_ene_int=0.0*splET2
cum_int=0.0
for i in range(size(splET2)-1):
    cum_int=cum_int+0.5*(splET2[i+1]+splET2[i])*(log_TK[i+1]-log_TK[i])
    spl_ene_int[i+1]=cum_int
    print (i)

OUTFILE=DIRPLOT+'spl_int_ene_over_t2'
np.savetxt(OUTFILE, np.c_[log_TK,spl_ene_int,splET2])

spl_free_energy=log_TK/log_TK[0]*EeV[0]-log_TK*spl_ene_int
spl_entropy=(splE-spl_free_energy)/log_TK

OUTFILE=DIRPLOT+'spl_free_ene_VS_temp'
np.savetxt(OUTFILE, np.c_[log_TK,spl_free_energy,spl_entropy])





#scipy.integrate.trapz(y, x=None, dx=1.0, axis=-1)



#scipy.integrate.trapz(ENE_)
# ax1.loglog(x,y,'-',color='r',linewidth=1.5,markersize=5,mew=1.5,mec='r',mfc='1.0',alpha=1.0,clip_on=True,label=r'$\alpha=1.0$')

#ax1.legend(loc=0,fontsize=18)
#plt.show()
#
#
#
#ax1.set_title(r'indbald dynamics $U_i/W=-0.5~U_f/W=-10.0 $')


# if args.outPDF:
#     figNAME=args.figNAME
#     figNAME=figNAME+'.pdf'
#     fig.savefig(figNAME,bbox_inches='tight')


# if args.outSVG:
#     figNAME=args.figNAME
#     figNAME=figNAME+'.svg'
#     fig.savefig(figNAME,bbox_inches='tight')





