from netCDF4 import Dataset
from astropy.io import fits
#from multiprocessing import Process
from tqdm import tqdm
from subprocess import check_call
import shutil, os

def single_walker(wid, data, msk):
  nstep, _ = msk.shape
  for i in tqdm(range(nstep)):
    # i1 records the most recent new state
    if msk[i,wid] == 1:
      i1 = i
      continue

    # copy previous state in chaine
    for key in data.variables.keys():
      if data[key].dimensions == ('time', 'x1', 'x2', 'x3'):
        data[key][i,:,:,wid] = data[key][i1,:,:,wid]
      elif data[key].dimensions == ('time', 'x2', 'x3'):
        data[key][i,:,wid] = data[key][i1,:,wid]
      elif data[key].dimensions == ('time', 'x3'):
        data[key][i,wid] = data[key][i1,wid]
      else:
        pass

def main_to_mcmc(fname, field = 'main'):
  #fname = 'vla_ideal_saturn-n1000'
  #os.remove('%s-mcmc.nc' % fname)
  #shutil.copy('%s-main.nc' % fname, '%s-tmp.nc' % fname)

  if field == 'main':
    data = Dataset('%s-main.nc' % fname, 'r+')
  else:
    data = Dataset('%s.%s.nc' % (fname, field), 'r+')
  msk = fits.open('%s.fits' % fname)[3].data
  nstep, nwalker = msk.shape

  data['time'][:] = range(nstep)
  data['time'].long_name = "steps"
  data['time'].units = "1"
  data['x2'][:] = range(len(data['x2'][:]))
  data['x2'].long_name = "model"
  data['x2'].units = "1"
  data['x3'][:] = range(len(data['x3'][:]))
  data['x3'].long_name = "walker"
  data['x3'].units = "1"

  pid = []
  for j in range(nwalker):
    print('- processing walker %d/%d' % (j,nwalker))
    single_walker(j, data, msk)

  data.close()

  if field == 'main':
    check_call('ncks -d time,0,%d %s-main.nc %s-mcmc.nc' %
      (nstep-1, fname, fname), shell = True)
    os.remove('%s-main.nc' % fname)
  else:
    check_call('ncks -d time,0,%d %s.%s.nc %s.%s-mcmc.nc' %
      (nstep-1, fname, field, fname, field), shell = True)
    os.remove('%s.%s.nc' % (fname, field))
