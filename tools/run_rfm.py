#! /usr/bin/env python3
from pylab import *
from netCDF4 import *
from collections import OrderedDict
import os, sys, glob, getopt, argparse

# search for the closest hit file
def SearchHitFile(folder, wmin, wmax, version = 'head'):
  hitfile, width = '', 1.E9
  for file in glob.glob(folder + '%s.*.hit' % version):
    param = os.path.basename(file).split('.')[1]
    wmin, wmax = map(float, param.split('-'))
    if (wmin >= wmin) and (wmax <= wmax) and (wmax - wmin < width):
      hitfile = file
      width = wmax - wmin
  return hitfile

# create hit file
def CreateHitFile(inp_hitran, wmin, wmax, version = ''):
  inp_wave    = '%d %d' % (wmin, wmax)
  out_hit     = '%s.%d-%d.hit' % (version, 0.5*wmin, 1.5*wmax)
  head        = '20%s %d %d' % (version, 0.5*wmin, 1.5*wmax)
  inp_hit     = '%s.%d-%d.inp' % (version, 0.5*wmin, 1.5*wmax)
  log_hit     = '%s.%d-%d.log' % (version, 0.5*wmin, 1.5*wmax)
  # check whether hit file exists
  if os.path.exists(out_hit):
    return out_hit
  print('# Generating new hit file: %s ...' % out_hit)
  with open(inp_hit, 'w') as file:
    file.write(inp_hitran + '\n')
    file.write(inp_wave + '\n')
    file.write(out_hit + '\n')
    file.write(head + '\n')
  os.system(args['hitbin'] + ' < %s > %s' % (inp_hit, log_hit))
  hitfile = out_hit
  return hitfile

# create rfm driver file
def CreateRfmDrv(driver, output):
  print('# Creating rfm.drv ...')
  with open(output, 'w') as file:
    for sec in driver:
      if driver[sec] != None:
        file.write(sec + '\n')
        file.write(' ' * 4 + driver[sec] + '\n')

# create rfm atmosphere file
def CreateRfmAtm(molecules, atm, output):
  print('# Creating rfm.atm ...')
  n_levels = atm.shape[0]
  n_molecules = len(molecules.split())
  with open(output, 'w') as file:
    file.write('%d\n' % n_levels)
    file.write('*HGT [km]\n')
    for i in range(n_levels):
      file.write('%.8g ' % atm['HGT'][i])
    file.write('\n*PRE [mb]\n')
    for i in range(n_levels):
      file.write('%.8g ' % atm['PRE'][i])
    file.write('\n*TEM [K]\n')
    for i in range(n_levels):
      file.write('%.8g ' % atm['TEM'][i])
    for i in range(n_molecules):
      name = molecules.split()[i]
      file.write('\n*' + name + ' [ppmv]\n')
      for j in range(n_levels):
        file.write('%.8g ' % atm[name][j])
    file.write('\n*END')

# run rfm model to get absorption coefficient table
def RunRfm(hitfile,
        wmin, wmax, wdel,
        tmin, tmax, tnum,
        molecules, atm_data, rundir = '.'):
  pwd = os.getcwd()
  # create rfm run file
  driver = OrderedDict([
      ('*HDR'  ,   'Header for rfm'),
      ('*FLG'  ,   'TAB CTM'),
      ('*SPC'  ,   '%.4f %.4f %.4f' % (wmin, wmax, wdel)),
      ('*GAS'  ,   molecules),
      ('*ATM'  ,   'rfm.atm'),
      ('*DIM'  ,   'PLV \n    %d %.4f %.4f' % (tnum, tmin, tmax)),
      ('*TAB'  ,   'tab_*.txt'),
      ('*HIT'  ,   os.path.join(pwd,hitfile)),
      ('*END'  ,   '')
      ])
  if not os.path.exists(rundir):
    os.makedirs(rundir)
  CreateRfmDrv(driver, rundir + '/rfm.drv')
  CreateRfmAtm(molecules, atm_data, rundir + '/rfm.atm')

  # run rfm
  print('# Running RFM in dir "%s"' % rundir)
  if rundir != '.':
    os.chdir(os.path.join(pwd,rundir))
    os.system('echo | %s > %s' % (os.path.join(pwd,args['rfm']), 'rfm.log'))
    os.chdir(pwd)
  else:
    os.system('echo | %s > %s' % (args['rfm'], 'rfm.log'))
  print('# RFM finished.')

# create input file to get absorption coefficient table
def CreateKcoeffInp(file,
        wmin, wmax, wdel,
        tmin, tmax, tnum,
        molecules, atm_data, rundir = '.'):
  file.write('# Molecular absorber\n')
  file.write('%d\n' % len(molecules.split()))
  file.write(molecules + '\n')
  file.write('# Molecule data files\n')
  for ab in molecules.split():
    file.write('%-40s\n' % (rundir + '/tab_'+ab.lower()+'.txt',))
  file.write('# Wavenumber range\n')
  file.write('%-14.6g%-14.6g%-14.6g\n' % (wmin,wmax,int((wmax - wmin)/wdel) + 1))
  file.write('# Relative temperature range\n')
  file.write('%-14.6g%-14.6g%-14.6g\n' % (tmin, tmax, tnum))
  file.write('# Number of vertical levels\n')
  file.write('%d\n' % len(atm_data['TEM']))
  file.write('# Temperature\n')
  for i in range(len(atm_data['TEM'])):
    file.write('%-14.6g' % atm_data['TEM'][-(i + 1)])
    if (i + 1) % 10 == 0:
      file.write('\n')
  if (i + 1) % 10 != 0:
    file.write('\n')
  file.write('# Pressure\n')
  for i in range(len(atm_data['PRE'])):
    file.write('%-14.6g' % atm_data['PRE'][-(i + 1)])
    if (i + 1) % 10 == 0:
      file.write('\n')
  if (i + 1) % 10 != 0:
    file.write('\n')

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('--hitbin',
      default = 'hitbin',
      help = 'hitbin exe'
      )
  parser.add_argument('--rfm',
      default = 'rfm',
      help = 'rfm exe'
      )
  parser.add_argument('--par',
      default = '',
      help = 'hitran par file, must specify if no hit'
      )
  parser.add_argument('--hit',
      default = '',
      help = 'hitran hit file, must specify if no par'
      )
  parser.add_argument('--atm',
      default = 'atm.txt',
      help = 'atmosphere profile'
      )
  parser.add_argument('--wave',
      default = '0. 0. 1',
      help = 'wavenumber grid'
      )
  parser.add_argument('--temp',
      default = '-5. 5. 3',
      help = 'temperature grid',
      )
  parser.add_argument('--molecule',
      default = '',
      help = 'molecular absorber'
      )
  parser.add_argument('--output',
      default = 'stdout',
      help = 'kcoeff output file'
      )
  parser.add_argument('--version',
      default = '12',
      help = 'Hitran version'
      )
  parser.add_argument('--rundir',
      default = '.',
      help = 'direcotry contains tab files'
      )

  # parse command-line inputs
  args = vars(parser.parse_args())

  # must specify either par or hit
  if args['par'] == '' and args['hit'] == '':
    raise Exception('Must specify either par or hit.')

  # atmospheric profile
  atm_data = genfromtxt(args['atm'], names = True)

  # Generate HITRON absorption coefficients
  wmin, wmax, wdel = map(float, args['wave'].split())
  tmin, tmax, tnum = map(float, args['temp'].split())
  tnum = int(tnum)
  if args['hit'] == '':
    hitfile = CreateHitFile(args['par'], wmin, wmax, version = args['version'])
  else:
    hitfile = args['hit']
  print("# Hit file is: ", hitfile)
  RunRfm(hitfile, wmin, wmax, wdel, tmin, tmax, tnum, args['molecule'], atm_data,
         rundir = args['rundir'])

  # Generate absorption coefficients input file
  if args['output'] == 'stdout':
    file = sys.stdout
  else:
    file = open(args['output'], 'w')
  CreateKcoeffInp(file, wmin, wmax, wdel, tmin, tmax, tnum, args['molecule'], atm_data,
                  rundir = args['rundir'])
  if file != sys.stdout:
    file.close()
