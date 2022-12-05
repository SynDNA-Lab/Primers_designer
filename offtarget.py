# NOTE: This file is a fucking mess, clean up the code
import os
import shutil
from itertools import groupby
from collections import Counter

import target
import bowtie
import config
from primer import Candidates



def generate_report():
    reppath = f"{config.homepath}/input/reports/{config.jobname}"
    if not os.path.isdir(reppath):
        os.mkdir(reppath)

    shutil.move("rmvduplic.csv", f"{reppath}/rmvduplic.csv")
    shutil.move("p3_result.txt", f"{reppath}/p3_result.txt")
    shutil.move("rmvchr.csv", f"{reppath}/rmvchr.csv")
    shutil.move("primercheck.fasta", f"{reppath}/primercheck.fasta")
    shutil.move("primercheck1.fasta", f"{reppath}/primercheck1.fasta")
    shutil.move("settings", f"{reppath}/settings")
    shutil.move("final_out.csv", f"{reppath}/final_out.csv")
    shutil.move("report.txt", f"{reppath}/report.txt")


def delete_multiple_element(list_object, indices):
    idc = sorted(list(set(indices)), reverse=True)
    for idx in idc:
        if idx < len(list_object):
            list_object.pop(idx)


def get_min_rid():
    records = target.load("primercheck1.fasta")
    ids = [r.id for r in records]

    global run_id_dict
    run_id_dict = {}
    for k, grp in groupby(ids, lambda x: x.split("_")[-3]):
        grp_list = list(grp)
        runs = [gl.split("_")[-2] for gl in grp_list]
        run_id_dict[k] = runs


def filter():
    # We will need to automatically create a bowtie index in the homepath/bowtie/index folder
    candall = [str(cand.target) for cand in Candidates.all ] # getting all candiate names 

    bowtie.run(index=config.jobname, inp_file="primercheck.fasta", out_file="rmvduplic.csv")

    btdata = bowtie.parse("rmvduplic.csv")

    bt_names = [bt[0] for bt in btdata if bt[0]]
    c = Counter([x.split("_")[-2] for x in bt_names])
    offtarget_idx =  [int(x) for x in list(Counter(el for el in c.elements() if c[el] > 3).keys())]
    delete_multiple_element(Candidates.all, offtarget_idx)
    bowtie.write_primers()

    bowtie.run(index="s_cerevisiae", inp_file="primercheck1.fasta", out_file="rmvchr.csv")
    btdata = bowtie.parse("rmvchr.csv")
    otnew = [x for x in btdata if x[2] != config.chromosome] # delte all off-targets that are on the same chromosome 
    # that was already scanned beforehand
    get_min_rid()
    sel_ids = chk_targ(otnew)
    
    btnewids = set([x[0].split("_")[-3] for x in otnew])
    btskipped = set(run_id_dict.keys()).difference(btnewids)
    sel_ids.extend([int(min(run_id_dict[ri])) for ri in list(btskipped)])

    sel_cand = [Candidates.all[s] for s in sel_ids]
    not_found = set(candall).difference(set([c.target for c in sel_cand]))

    with open("final_out.csv", "w") as file:
        file.write("Target,Pair_Penalty,FWD,REV,FWD_Length,REV_Length,FWD_TM,REV_TM,FWD_GC,REV_GC,PCR_Size\n")
        for sel in sel_cand:
            file.write(f"{sel.target},{sel.pair_penalty},{sel.left_sequence},{sel.right_sequence},{len(sel.left_sequence)},{len(sel.right_sequence)},{sel.left_tm},{sel.right_tm},{sel.left_gc_percent},{sel.right_gc_percent},{sel.pair_product_size}\n")

    with open(f"report.txt", "w") as file:
        file.write("Failed to find primers for\n")
        file.write(str(not_found))


def has_offtarget(l1):
    if len(l1) < 2:
        return False
    fwd_facing = [int(l[3]) for l in l1 if l[1] == "+"]
    rev_facing = [int(l[3]) for l in l1 if l[1] == "-"]
    for rev in rev_facing:
        for fwd in fwd_facing:
            if (rev - fwd > 0) and (rev - fwd <= 10_000):
                return True # Offtarget w/ less then 10kb found
    return False # No offtarget


def chk_chr(pid_list):
    for _, chr_group in groupby(pid_list, lambda x: x[2]):
        ot_chr = list(chr_group)
        if has_offtarget(ot_chr):
            return False
    return True


def chk_id(trg_list, clist):
    for k1, pid_group in groupby(trg_list, lambda x: x[0].split("_")[-2]):
        pid_list = list(pid_group)
        current_id = pid_list[0][0].split("_")[-2]

        if clist:
            if current_id > min(clist):
                return True, min(clist)
        if len(pid_list) < 2:
            return True, current_id
        elif len(pid_list) > 5:
            continue
        if chk_chr(pid_list):
            return True, current_id
    return False, None


def chk_targ(offtargets):
    b1 = True
    sel = []
    for key, trg_group in groupby(offtargets, lambda x: x[0].split("_")[-3]):
        trg_list = list(trg_group)
        if b1:
            b1 = False
            print(trg_list)

        prid = [t[0].split("_")[-2] for t in trg_list]
        chckidbuffer = set(run_id_dict[key]).difference(set(prid))

        x, y = chk_id(trg_list, chckidbuffer)
        if x:
            sel.append(int(y))
        else:
            print("-------")
            print(trg_list[0][0]) # Error, no suitable primer found
            print("-------")

    return sel