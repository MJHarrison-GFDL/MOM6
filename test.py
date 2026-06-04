import numpy as np
import regrid
#import remap
import matplotlib.pyplot as plt#import matplotlib.pyplot as plt


def write_MOM_PF(ale_coord_config="UNIFORM",ale_units="m",eos='LINEAR',interpolation_scheme='PLM',p_ref=0.e7):

    file=open('MOM_input','w')
    txt="REGRIDDING_COORDINATE_UNITS = \""+ale_units+" \"\n"
    print(txt)
    file.write(txt)
    txt="ALE_COORDINATE_CONFIG = \""+ale_coord_config+" \"\n"
    file.write(txt)
    txt="EQN_OF_STATE = \""+eos+" \"\n"
    file.write(txt)
    txt="INTERPOLATION_SCHEME = \""+interpolation_scheme+" \"\n"
    file.write(txt)
    txt="P_REF = "+str(p_ref)+" \n"
    file.write(txt)
    file.close()

def destroy_MOM_PF():
    import os
    if os.path.isfile('MOM_input'):
        os.remove('MOM_input')
    if os.path.isfile('MOM_parameter_doc.short'):
        os.remove('MOM_parameter_doc.short')
    if os.path.isfile('MOM_parameter_doc.all'):
        os.remove('MOM_parameter_doc.all')
    if os.path.isfile('MOM_parameter_doc.layout'):
        os.remove('MOM_parameter_doc.layout')
    if os.path.isfile('MOM_parameter_doc.debugging'):
        os.remove('MOM_parameter_doc.debugging')
def wright_eos(T,S,p):
  """

 **********************************************************************
   The subroutines in this file implement the equation of state for   *
   sea water using the formulae given by  Wright, 1997, J. Atmos.     *
   Ocean. Tech., 14, 735-740.  Coded by R. Hallberg, 7/00.            *
 ***********************************************************************
    Converted to Python from F90 by M Harrison 10/11.
 Calculate seawater equation of state, given T[degC],S[PSU],p[Pa]
 Returns density [kg m-3]

 ***********************************************************************

 """

  a0 = 7.057924e-4; a1 = 3.480336e-7; a2 = -1.112733e-7;
  b0 = 5.790749e8;  b1 = 3.516535e6;  b2 = -4.002714e4;
  b3 = 2.084372e2;  b4 = 5.944068e5;  b5 = -9.643486e3;
  c0 = 1.704853e5;  c1 = 7.904722e2;  c2 = -7.984422;
  c3 = 5.140652e-2; c4 = -2.302158e2; c5 = -3.079464;

  al0 = a0 + a1*T +a2*S
  p0  = b0 + b4*S + T * (b1 + T*(b2 + b3*T) + b5*S)
  lam = c0 +c4*S + T * (c1 + T*(c2 + c3*T) + c5*S)
  I_denom = 1.0 / (lam + al0*(p+p0))
  rho = (p + p0) * I_denom

  return rho

VERBOSE=False
H=100.;n0=10;z0=np.linspace(0,-H,n0+1)
T_s=2;T_b=-2
S_s=34.2;S_b=35.
T=np.linspace(T_s,T_b,n0)
S=np.linspace(S_s,S_b,n0)
ps=np.array([1.e4])
fs=np.array([0.])
zbot=np.array([H])
nfig=1

print('====== ZSTAR (UNIFORM) ======')
#Uniform Resolution
CoordVals=z0
#Perturb the Surface
z0[0]=z0[0]+0.1
z_in = z0[np.newaxis,np.newaxis,:]
T_in = T[np.newaxis,np.newaxis,:]
S_in = S[np.newaxis,np.newaxis,:]
z_out = np.zeros(z_in.shape)
ps_in = ps[np.newaxis,:]
fs_in = fs[np.newaxis,:]
zbot_in=zbot[np.newaxis,:]
destroy_MOM_PF()
write_MOM_PF()
z_out=regrid.regrid_mod.update_grid(z_in,T_in,S_in,zbot,ps_in,fs_in,'Z*',CoordVals,'PCM',False)
diff=z_out-z_in
print('RMS motion= ',np.std(diff))
print('Initial Interface Positions=',np.squeeze(z_in))
print('Final Interface Positions=',np.squeeze(z_out))
if VERBOSE:
    plt.figure(nfig);nfig=nfig+1
    plt.plot(np.squeeze(z_in).T,'bo')
    plt.plot(np.squeeze(z_out).T,'rx')
    plt.title('Interface Positions: Initial(o);Final(x)')

print('====== ZSTAR (FNC1) ======')
CoordVals=z0 # This is ignored in the subsequent call to update_grid.
#Perturb the Surface
z0[0]=z0[0]+0.1
z_in = z0[np.newaxis,np.newaxis,:]
T_in = T[np.newaxis,np.newaxis,:]
S_in = S[np.newaxis,np.newaxis,:]
z_out = np.zeros(z_in.shape)
ps_in = ps[np.newaxis,:]
fs_in = fs[np.newaxis,:]
zbot_in=zbot[np.newaxis,:]
destroy_MOM_PF()
write_MOM_PF(ale_coord_config="FNC1:1.0,100.,0,1.e-6",ale_units="m")
z_out=regrid.regrid_mod.update_grid(z_in,T_in,S_in,zbot,ps_in,fs_in,'Z*',CoordVals,'PCM',False)
print('Initial Interface Positions=',np.squeeze(z_in))
print('Final Interface Positions=',np.squeeze(z_out))
diff=z_out-z_in
print('RMS motion= ',np.std(diff))


if VERBOSE:
    plt.figure(nfig);nfig=nfig+1
    plt.plot(np.squeeze(z_in).T,'bo')
    plt.plot(np.squeeze(z_out).T,'rx')
    plt.title('Interface Positions: Initial(o);Final(x)')


print('====== RHO (Linear EOS;RFNC1) ======')
RHO_T0_S0=1000.
dRdT=-0.2
dRdS=0.8
rho=RHO_T0_S0+dRdT*T + dRdS*S
CoordVals=np.linspace(rho[0],rho[-1],n0+1)
print(CoordVals)
#Perturb the Surface
z0[0]=z0[0]+0.1
z_in = z0[np.newaxis,np.newaxis,:]
T_in = T[np.newaxis,np.newaxis,:]
S_in = S[np.newaxis,np.newaxis,:]
z_out = np.zeros(z_in.shape)
ps_in = ps[np.newaxis,:]
fs_in = fs[np.newaxis,:]
zbot_in=zbot[np.newaxis,:]
destroy_MOM_PF()
write_MOM_PF(ale_coord_config="\"RFNC1:10,1026.96,1027.12,1027.28,1.065,1028.4,0.001,1\"",ale_units="kg m^-3",eos='LINEAR',interpolation_scheme='PLM')
#write_MOM_PF(ale_coord_config="\"RFNC1:10,1026.96,1027.12,1027.28,1.12,1028.4,0.001,0\"",ale_units="kg m^-3",eos='LINEAR',interpolation_scheme='PLM')
z_out=regrid.regrid_mod.update_grid(z_in,T_in,S_in,zbot,ps_in,fs_in,'RHO',CoordVals,'PCM',True)
diff=z_out-z_in
print('Initial Interface Positions=',np.squeeze(z_in))
print('Final Interface Positions=',np.squeeze(z_out))
print('RMS motion= ',np.std(diff))

if VERBOSE:
#    print('Initial Interface Positions=',np.squeeze(z_in))
#    print('Final Interface Positions=',np.squeeze(z_out))
    plt.figure(nfig);nfig=nfig+1
    plt.plot(np.squeeze(z_in).T,'bo')
    plt.plot(np.squeeze(z_out).T,'rx')
    plt.title('Interface Positions: Initial(o);Final(x)')


print('====== RHO (WRIGHT;RFNC1) ======')
rho=wright_eos(T,S,p=2.e7)
#Perturb the Surface
z0[0]=z0[0]+0.1
z_in = z0[np.newaxis,np.newaxis,:]
T_in = T[np.newaxis,np.newaxis,:]
S_in = S[np.newaxis,np.newaxis,:]
z_out = np.zeros(z_in.shape)
ps_in = ps[np.newaxis,:]
fs_in = fs[np.newaxis,:]
zbot_in=zbot[np.newaxis,:]
destroy_MOM_PF()
write_MOM_PF(ale_coord_config="RFNC1:10,1036.5,1036.67,1036.8,0.8,1037.62,0.1,1",ale_units="kg m^-3",eos='WRIGHT',interpolation_scheme='PLM',p_ref=2.e7)
z_out=regrid.regrid_mod.update_grid(z_in,T_in,S_in,zbot,ps_in,fs_in,'RHO',CoordVals,'PCM',False)
diff=z_out-z_in
print('RMS motion= ',np.std(diff))
print('Initial Interface Positions=',np.squeeze(z_in))
print('Final Interface Positions=',np.squeeze(z_out))

if VERBOSE:
#    print('Initial Interface Positions=',np.squeeze(z_in))
#    print('Final Interface Positions=',np.squeeze(z_out))
    plt.figure(nfig);nfig=nfig+1
    plt.plot(np.squeeze(z_in).T,'bo')
    plt.plot(np.squeeze(z_out).T,'rx')
    plt.title('Interface Positions: Initial(o);Final(x)')

# print('====== HYCOM1 (WRIGHT;RFNC1) ======')
# rho=wright_eos(T,S,p=2.e7)
# CoordRes=rho
# print(CoordRes)
# #Perturb the Surface
# z0[0]=z0[0]+0.1
# z_in = z0[np.newaxis,np.newaxis,:]
# T_in = T[np.newaxis,np.newaxis,:]
# S_in = S[np.newaxis,np.newaxis,:]
# z_out = np.zeros(z_in.shape)
# ps_in = ps[np.newaxis,:]
# fs_in = fs[np.newaxis,:]
# zbot_in=zbot[np.newaxis,:]
# destroy_MOM_PF()
# write_MOM_PF(ale_coord_config="HYBRID:hycom1_10.nc,sigma2",ale_units="kg m^-3",eos='WRIGHT',interpolation_scheme='PLM',p_ref=2.e7)
# z_out=regrid.regrid_mod.update_grid(z_in,T_in,S_in,zbot,ps_in,fs_in,'HYCOM1',CoordRes,'PCM')
# diff=z_out-z_in
# print('RMS motion= ',np.std(diff))

# if VERBOSE:
# #    print('Initial Interface Positions=',np.squeeze(z_in))
# #    print('Final Interface Positions=',np.squeeze(z_out))
#     plt.figure(nfig);nfig=nfig+1
#     plt.plot(np.squeeze(z_in).T,'bo')
#     plt.plot(np.squeeze(z_out).T,'rx')
#     plt.title('Interface Positions: Initial(o);Final(x)')



#destroy_MOM_PF()
if VERBOSE:
    plt.show()
