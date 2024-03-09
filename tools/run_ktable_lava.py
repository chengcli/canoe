#! python3
from multiprocessing import Pool, cpu_count
import os, re, subprocess

cmake_source_dir = "/home/linfel/canoe"
cmake_binary_dir = "/home/linfel/canoe/build"

# file architecture
hitbin = f"{cmake_binary_dir}/bin/hitbin.release"
rfm = f"{cmake_binary_dir}/bin/rfm.release"
run_rfm = f"{cmake_binary_dir}/bin/run_rfm.py"
kcoeff = f"{cmake_binary_dir}/bin/kcoeff.release"
cktable = f"{cmake_binary_dir}/bin/run_cktable.py"

hitran = ""
# use .par file explicitly as hitfile
hitfile = "/home/linfel/co_HITRAN2020/lb_28Si-16O__SiOUVenIR__100-110__296K.par"

# number of threads
max_threads = cpu_count() // 2

# ktable specifics
atm = f"{cmake_source_dir}/data/lava_SiO_atm_isothermal.txt"
# mol     = "CH4 C2H2 C2H4 C2H6"
mol = "SiO"

# temperature grid
temp = "-5 5 3"

# flags
generate_tab = True
generate_nc = True
generate_cktable = False

# continuum absorption
cia = "H2-H2-eq H2-He-eq"

# spectral bands
band = []
band.append("100 110 0.01")
# band.append("200 400 0.01")
# band.append("400 600 0.01")
# band.append("600 700 0.01")
# band.append("700 860 0.01")
# band.append("860 960 0.01")
# band.append("960 1200 0.01")
# band.append("1200 1400 0.01")
# band.append("1400 1700 0.01")
# band.append("1700 2100 0.01")
# band.append("2100 3450 0.01")
# band.append("3450 4800 0.01")
# band.append("4800 6300 0.01")
# band.append("6300 7800 0.01")
# band.append("7800 9300 0.01")

# number of parallel threads
nthreads = min(max_threads, len(band))


# run ktable in single thread
def RunSingleKtable(wave):
    wname = re.sub(" ", "-", wave)
    print("working on band %s ..." % wname)
    wmin, wmax, dw = map(float, wave.split())
    tab_folder = str(wmin) + "-" + str(wmax)
    kinp = "kcoeff.inp-" + wname
    kncfile = "kcoeff." + wname + ".nc"
    # create tab files and kcoeff.inp
    if generate_tab:
        script = [
            "python",
            run_rfm,
            "--hitbin",
            hitbin,
            "--rfm",
            rfm,
            "--par",
            hitran,
            "--hit",
            hitfile,
            "--atm",
            atm,
            "--wave",
            wave,
            "--temp",
            temp,
            "--molecule",
            mol,
            "--output",
            kinp,
            "--rundir",
            tab_folder,
        ]
        out, err = subprocess.Popen(
            script, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        ).communicate()
        print(out.decode(), err.decode())
        # if not ('Successful' in out.decode()):
        #  raise RuntimeError("Error in generating tab files.")

    # run kcoeff
    if generate_nc:
        script = [kcoeff, "-i", kinp, "-o", "kcoeff.%s.nc" % wname]
        out, err = subprocess.Popen(
            script, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        ).communicate()
        print(out.decode(), err.decode())
        # if err != b'':
        #    raise RuntimeError("Error in generating kcoeff.**.nc.")

    # create correlated-K table
    if generate_cktable:
        script = [
            cktable,
            "--kcoeff",
            kncfile,
            "--atm",
            atm,
            "--cia",
            cia,
            "--input",
            inpfile,
        ]
        out, err = subprocess.Popen(
            script, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        ).communicate()
        print(out.decode(), err.decode())
    print("band %s finishes." % wave)


if __name__ == "__main__":
    # parallel on spectral bands
    pool = Pool(nthreads)
    pool.map(RunSingleKtable, band)
