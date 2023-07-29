#! /usr/bin/env python3
import argparse, os, shutil
from copy import copy
from glob import glob
from datetime import datetime
from subprocess import check_call
from numpy import sort
from main2mcmc import main_to_mcmc


def CombineTimeseries(case, field, stamps, path="./", remove=False):
    print("Concatenating output field %s ..." % field)
    fname = ""
    for stamp in stamps:
        fname += path + "/%s.%s.%s.nc " % (case, field, stamp)
    target = "%s.%s.nc" % (case, field)

    if len(stamps) > 1:
        # check_call('ncrcat -h %s -o %s' % (fname, target), shell = True)
        check_call(
            "ncrcat -h %s/%s.%s.?????.nc -o %s/%s" % (path, case, field, path, target),
            shell=True,
        )
        if remove:
            for f in fname.split():
                os.remove(f)
    else:
        shutil.move(fname[:-1], target)


def CombineFields(case, fields, output, path="./"):
    # if more than two fields, a combine_rules file is needed
    flist = sort(list(set([f[3:] for f in fields])))

    if not os.path.isfile("combine_rules"):
        os.system('echo "%s -> main" > combine_rules' % ",".join(flist))
    with open("combine_rules", "r") as file:
        for line in file.readlines():
            fid = list(map(int, line.split()[0].split(",")))
            name = line.split()[2]

            print("Combining output fields: ", fid, "to", name)
            if output != "None":
                ncout = "%s/%s-%s-%s.nc" % (path, case, output, name)
            else:
                ncout = "%s/%s-%s.nc" % (path, case, name)

            if len(fid) > 1:  # combining
                fname1 = path + "/%s.out%d.nc" % (case, fid[0])
                for i in fid[1:]:
                    fname2 = path + "/%s.out%d.nc" % (case, i)
                    check_call("ncks -A %s %s" % (fname2, fname1), shell=True)
                    os.remove(fname2)
            else:  # renaming
                fname1 = path + "/%s.out%d.nc" % (case, fid[0])
            shutil.move(fname1, ncout)


def ParseOutputFields(path):
    files = glob(path + "/*.out*.[0-9][0-9][0-9][0-9][0-9].nc")
    cases = []
    fields = []
    stamps = []

    for fname in files:
        field, stamp, ext = os.path.basename(fname).split(".")[-3:]
        case = ".".join(os.path.basename(fname).split(".")[:-3])
        if case not in cases:
            cases.append(case)
        if field not in fields:
            fields.append(field)
        if stamp not in stamps:
            stamps.append(stamp)
    stamps = sorted(stamps)
    return cases, fields, stamps


def CombineFITS(case, output, path="./", remove=False):
    print("Combining FITS output ...")
    files = glob(path + "/%s.out?.[0-9][0-9][0-9][0-9][0-9].fits" % case)
    if len(files) == 0:
        return None
    if output != "None":
        fitsout = "%s-%s.fits" % (case, output)
    else:
        fitsout = "%s.fits" % case
    root = os.path.dirname(os.path.realpath(__file__))
    check_call(
        "%s/bin/fitsmerge -o %s -i %s/%s.out?.?????.fits" % (root, fitsout, path, case),
        shell=True,
    )
    if remove:
        for f in files:
            os.remove(f)
    return fitsout


if __name__ == "__main__":
    now = datetime.now()
    today = "%02d%02d" % (now.month, now.day)

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d", "--dir", default=".", help="directory of the simulation to combine"
    )
    parser.add_argument(
        "-o", "--output", default="None", help="appending additional name to the output"
    )
    parser.add_argument(
        "--no-remove", action="store_true", help="do not remove original files"
    )
    parser.add_argument(
        "--no-merge", action="store_true", help="do not merge different fields"
    )
    parser.add_argument(
        "--no-main2mcmc", action="store_true", help="do not perform main2mcmc change"
    )
    args = vars(parser.parse_args())

    cases, fields, stamps = ParseOutputFields(args["dir"])
    fields.sort()

    print("##########################")
    print("## Combine output files ##")
    print("##########################")
    for case in cases:
        print("Working on case %s..." % case)
        fitsout = CombineFITS(
            case, args["output"], path=args["dir"], remove=not args["no_remove"]
        )
        for field in fields:
            CombineTimeseries(
                case, field, stamps, remove=not args["no_remove"], path=args["dir"]
            )
        if not args["no_merge"]:
            CombineFields(case, fields, args["output"], path=args["dir"])
        if not args["no_main2mcmc"] and fitsout:
            if args["no_merge"]:
                for field in fields:
                    main_to_mcmc(fitsout[:-5], fields)
            else:
                main_to_mcmc(fitsout[:-5])
    print("Done.\n")
