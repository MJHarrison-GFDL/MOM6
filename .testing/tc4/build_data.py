import netCDF4 as nc
import numpy as np

x = nc.Dataset('ocean_hgrid.nc').variables['x'][1::2, 1::2]
y = nc.Dataset('ocean_hgrid.nc').variables['y'][1::2, 1::2]
zbot = nc.Dataset('topog.nc').variables['depth'][:]
zbot0 = zbot.max()


def t_fc(x, y, z, radius=5.0, tmag=1.0):
    """a radially symmetric anomaly in the center of the domain.
    units are meters and degC.
    """
    ny, nx = x.shape
    nz = z.shape[0]

    x0 = x[int(ny/2), int(nx/2)]
    y0 = y[int(ny/2), int(nx/2)]

    tl = np.zeros((nz, ny, nx))
    zb = z[-1]
    if len(z) > 1:
        zd = z / zb
    else:
        zd = [0.]
    for k in np.arange(len(zd)):
        r = np.sqrt((x - x0)**2 + (y - y0)**2)
        tl[k, :] += (1.0 - np.minimum(r / radius, 1.0)) * tmag * (1.0 - zd[k])
    return tl



### Create an initial condition/sponge file . This file contains interface
### positions for the target data and the cell average values for temperature
### and salinity

ny, nx = x.shape
nz = 3
z = (np.arange(nz) * zbot0) / nz
zi = np.zeros(nz+1)
zi[0]=0.
for k in np.arange(1,nz):
    zi[k]=0.5*(z[k]+z[k-1])
zi[-1]=zbot0
Eta=zi[np.newaxis,:,np.newaxis,np.newaxis]
Eta=np.tile(Eta,(1,1,ny,nx))

temp = t_fc(x, y, z)
salt = np.zeros(temp.shape)+35.0
fl = nc.Dataset('temp_salt_ic.nc', 'w', format='NETCDF3_CLASSIC')
fl.createDimension('lon', nx)
fl.createDimension('lat', ny)
fl.createDimension('z', nz)
fl.createDimension('zi', nz+1)
fl.createDimension('Time', None)
zv = fl.createVariable('z', 'f8', ('z'))
ziv = fl.createVariable('zi', 'f8', ('zi'))
lonv = fl.createVariable('lon', 'f8', ('lon'))
latv = fl.createVariable('lat', 'f8', ('lat'))
timev = fl.createVariable('Time', 'f8', ('Time'))
timev.calendar = 'noleap'
timev.units = 'days since 0001-01-01 00:00:00.0'
timev.modulo = ' '
tv = fl.createVariable('ptemp', 'f8', ('Time', 'z', 'lat', 'lon'),
                       fill_value=-1.e20)
sv = fl.createVariable('salt', 'f8', ('Time', 'z', 'lat', 'lon'),
                       fill_value=-1.e20)
etav = fl.createVariable('Eta', 'f8', ('Time', 'zi', 'lat', 'lon'),
                       fill_value=-1.e20)
tv[:] = temp[np.newaxis, :]
sv[:] = salt[np.newaxis, :]
etav[:] = Eta
zv[:] = z
ziv[:] = zi
lonv[:] = x[0, :]
latv[:] = y[:, 0]
timev[0] = 0.
fl.sync()
fl.close()


### Construct the target values for u/v which are identically
### zero in this case.

xu = nc.Dataset('ocean_hgrid.nc').variables['x'][1::2, 0::2]
yu = nc.Dataset('ocean_hgrid.nc').variables['y'][1::2, 0::2]
xv = nc.Dataset('ocean_hgrid.nc').variables['x'][0::2, 1::2]
yv = nc.Dataset('ocean_hgrid.nc').variables['y'][0::2, 1::2]

(ny,nxp)=xu.shape
(nyp,nx)=xv.shape
u = np.zeros((nz,ny,nxp))
v = np.zeros((nz,nyp,nx))

fl = nc.Dataset('uv_sponge.nc', 'w', format='NETCDF3_CLASSIC')
fl.createDimension('lonu', nxp)
fl.createDimension('latu', ny)
fl.createDimension('lonv', nx)
fl.createDimension('latv', nyp)
fl.createDimension('z', nz)
fl.createDimension('Time', None)
zv = fl.createVariable('z', 'f8', ('z'))
lonuv = fl.createVariable('lonu', 'f8', ('lonu'))
latuv = fl.createVariable('latu', 'f8', ('latu'))
lonvv = fl.createVariable('lonv', 'f8', ('lonv'))
latvv = fl.createVariable('latv', 'f8', ('latv'))
timev = fl.createVariable('Time', 'f8', ('Time'))
timev.calendar = 'noleap'
timev.units = 'days since 0001-01-01 00:00:00.0'
timev.modulo = ' '
uv = fl.createVariable('UVEL', 'f8', ('Time', 'z', 'latu', 'lonu'),
                       fill_value=-1.e20)
vv = fl.createVariable('VVEL', 'f8', ('Time', 'z', 'latv', 'lonv'),
                       fill_value=-1.e20)
uv[:] = u[np.newaxis, :]
vv[:] = v[np.newaxis, :]
zv[:] = z
lonuv[:] = xu[0, :]
latuv[:] = yu[:, 0]
lonvv[:] = xv[0, :]
latvv[:] = yv[:, 0]
timev[0] = 0.
fl.sync()
fl.close()


# Make Sponge forcing file

dampTime = 20.0     # days
secDays = 8.64e4
fl = nc.Dataset('sponge.nc', 'w', format='NETCDF3_CLASSIC')
fl.createDimension('lon', nx)
fl.createDimension('lat', ny)
lonv = fl.createVariable('lon', 'f8', ('lon'))
latv = fl.createVariable('lat', 'f8', ('lat'))
spv = fl.createVariable('Idamp', 'f8', ('lat', 'lon'), fill_value=-1.e20)
Idamp = np.zeros((ny, nx))
if dampTime > 0.:
    Idamp = 0.0 + 1.0 / (dampTime * secDays)
spv[:] = Idamp
lonv[:] = x[0, :]
latv[:] = y[:, 0]
fl.sync()
fl.close()
