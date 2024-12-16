import imp # Used to import general.py
import time # Used for time
import numpy as np
import pandas as pd
import itertools as it
import subprocess
from subprocess import Popen

imp.load_source("general",
                ("/lab/bartel1_ata/mcgeary/computation/"
                 "AgoRBNS/general/general.py"
                )
               )
imp.load_source("RBNS_methods",
                ("/lab/bartel1_ata/mcgeary/computation/"
                 "AgoRBNS/general/RBNS_methods.py"
                )
               )
from general import *
from RBNS_methods import *
from sitetypes import get_seq_site_map

def bjobs():
    return(Popen("bjobs", shell=True,
          stdout=PIPE, stderr=PIPE).communicate()[0].split("\n"))[1:-1]

# def count_lines(path):
#     return(int(Popen("wc -l %s" %(path), shell=True,
#           stdout=PIPE).communicate()[0].split()[0]))

# def ls(path):
#     return(Popen("ls %s" %(path), shell=True,
#           stdout=PIPE, stderr=PIPE).communicate()[0].strip().split())

def check_jobs():
    return(ls("job_queue/*_jobs.txt"))

def check_master():
    return(ls("job_queue/*_master.txt"))


class Job:
    def __init__(self, job_id, status, host):
        self.id = job_id
        self.status = status
        self.host = host
        self.len = 1
    def __len__(self):
        return self.len
    def __add__(self, other):
        if other.__class__.__name__ == "str" and other == self.host:
            self.len += 1
            return(self)
    def __radd__(self, other):
        if other == 0:
            return(len(self))
        elif other.__class__.__name__ == "Job":
            return(len(self) + len(other))
    def __iadd__(self, other):
        if other.__class__.__name__ == "str" and other == self.host:
            self.len += 1
            return(self)
    def __repr__(self):
        if len(self) == 0:
            return("%s\t%s" %(self.id, self.status))
        else:
            return("%s: %s %s %s" %(self.id, self.status, self.host,
                                 len(self)))

class JobList():
    def __init__(self):
        unfinished = True
        while unfinished == True:
            job_ids = get_job_ids()
            hosts = get_job_hosts()
            status = get_job_statuses()
            if len(job_ids) == len(hosts) and len(job_ids) == len(status):
                check = [i_2 for i, j in enumerate(job_ids)
                              for i_2 in [hosts[i], status[i]]
                              if j != ""]
                if "" not in check:
                    unfinished = False
        job_list = []
        if len(job_ids) > 0:
            job_tuples = zip(job_ids, status, hosts)
            _job = Job(job_tuples[0][0], job_tuples[0][1], job_tuples[0][2])

            for job_id, status, host in job_tuples[1:]:
                # print(job_id)
                # print(status)
                # print(host)
                if job_id != "":
                    job_list.append(_job)
                    # print([job_id, status, host])
                    _job = Job(job_id, status, host)
                # print(_job)
                _job += host
            job_list.append(_job)
        self.jobs = job_list
        self.ids = [_job.id for _job in self.jobs]
        self.container = {}
        for _job in self.jobs:
            self.container[_job.id] = _job
    def __getitem__(self, key):
        if key in self.container:
            return(self.container[key])
        else:
            return
    def __len__(self):
        return(sum([len(_job) for _job in self.jobs]))



def get_all_hosts():
    print(" in get_all_hosts")
    hosts = [i.split()[0] for i in Popen("bhosts", shell=True,
             stdout=PIPE).communicate()[0].strip().split("\n")[1:]]
    print(" got hosts")
    return(hosts)

def get_job_info(name):
    titles = Popen("bhosts | head -1", shell=True,
                 stdout=PIPE).communicate()[0].strip().split()
    host_info = Popen("bhosts | grep %s" %(name), shell=True,
                 stdout=PIPE).communicate()[0].strip().split()
    return {key : value for key, value in zip(titles, host_info)[1:]}



class Host():
    def __init__(self, name):
        print("initializing Host object.")
        self.name = name
        print("Done initializing Host object.")
    def get_jobs(self):
        job_info = get_job_info(self.name)
        return(int(job_info["MAX"]) - int(job_info["RUN"]))
    def get_status(self):
        job_info = get_job_info(self.name)
        return(job_info["STATUS"])
    def __getitem__(self, key):
        if key == "jobs":
            return(self.get_jobs())
        elif key == "status":
            return(self.get_status())
        else:
            return ValueError

class HostList():
    EXCLUDED_HOSTS = ["it-c13b01", "it-c13b02", "it-c13b03", "it-c13b04",
                      "it-c13b05", "it-c13b06", "it-c13b07", "it-c13b08",
                      "it-c13b09", "it-c13b10", "it-c13b12",
                      "it-c13b13", "it-c13b14", "it-c13b15", "it-c13b16",
                      "it-c01b10", "it-c03b10", "it-c06b14", "it-r16u13",
                      "rou-c01b08", "rou-c05b09", "rou-c06b15", "rou-c06b16",
                      "airstream"]
    def __init__(self):
        print("Initializing HostList object")
        print("About to use 'get_al_hosts()'")
        self.list = get_all_hosts()
        print("Just got all hosts; about to assign self.hosts objects.")
        self.hosts = {name : Host(name) for name in self.list}
        print("Done assigning self.hosts objects; leaving initialization function.")
    def __len__(self):
        len_out = 0
        for name in self.list:
            print(name)
            print(self.hosts[name].get_status())
            print(self.hosts[name].get_jobs())
            print(name not in HostList.EXCLUDED_HOSTS)
            if (
                    self.hosts[name].get_status() == "ok" and
                    self.hosts[name].get_jobs() >= 0 and
                    name not in HostList.EXCLUDED_HOSTS
                ):
                len_out += 1
        return(len_out)
    def __getitem__(self, key):
        if key == "active":
            print("in active thing")
            names = self.list
            print("got host names in active")
            hosts = [(name, self.hosts[name]["jobs"],
                      self.hosts[name]["status"])
                     for name in names]
            print("got hosts")
            hosts = [host for host in hosts
                     if host[0] not in HostList.EXCLUDED_HOSTS
                     and host[1] >= 0 and host[2] == "ok"]
            print(hosts)
            return(sorted(hosts, key = lambda host: host[1], reverse=True))

def main():
    time_start = time.time()
    arguments = ["job", "-test_job_binary"]
    args = parse_arguments(arguments)
    jobs = args[0]
    test_job = args[1]
    if jobs == "hosts":
        hosts = get_all_hosts()
        return
    elif jobs == "restart":
        restart = True
    else:
        restart = False
    if restart:
        job_path = ls("job_queue/")
        script_id = job_path[0].split("_master")[0]
        master_queue_path = "job_queue/%s_master.txt" %(script_id)
        job_queue_path = "job_queue/%s_jobs.txt" %(script_id)
        with open(master_queue_path, "r+") as file_open:
            jobs = [i.strip() for i in list(file_open)]
    else:
        jobs = jobs.strip().split("\\n")
        if jobs[-1] == "":
            jobs = jobs[:-1]
        script_id = "%s_%s" %(time.time(), random.random())
        master_queue_path = "job_queue/%s_master.txt" %(script_id)
        job_queue_path = "job_queue/%s_jobs.txt" %(script_id)
    if (check_master() and not restart) or test_job:
        print("giving jobs")
        with open(job_queue_path, "w+") as queue_file:
            queue_file.write("\n".join(jobs))
        accepted = False
        while not accepted:
            time.sleep(1)
            completed = ls("job_queue/%s_read.txt" %(script_id))
            if completed:
                os.remove(job_queue_path)
                os.remove("job_queue/%s_read.txt" %(script_id))
                accepted = True
                return
    else:
        print("this is master")
        with open(master_queue_path, "w+") as master_file:
            master_file.write("\n".join(jobs))
        print("wrote master file")
        prior_host = ""
        host = ""
        jobs_accumulated = 0
        jobs_completed = []
        jobs_active = {}
        say_finished = True
        while jobs or jobs_active:
            outside_jobs = check_jobs()
            while outside_jobs:
                for path in outside_jobs:
                    with open(path, "r") as file_in:
                        jobs_new = file_in.read().split("\n")
                    jobs += jobs_new
                    read_path = path.split("_jobs.txt")[0] + "_read.txt"
                    with open(read_path, "w+") as read_out:
                        print("got jobs")
                time.sleep(5)
                outside_jobs = check_jobs()
            if len(jobs) > 0:
                print("updating master file")
                with open(master_queue_path, "w+") as master_file:
                    master_file.write("\n".join(jobs))
                job = jobs[0]
                mirna, exp, cond = job.split()[2:5]
                mirna = mirna.split(",")[0]
                if exp.split("_")[-1] == "nb":
                    nb = True
                    tp = False
                elif exp.split("_")[-1] == "tp":
                    nb = False
                    tp = True
                else:
                    nb = False
                    tp = False
                cond_split = cond.split("_")
                if len(cond_split) == 2 and exp != "equil_flowthrough":
                    cond, rep = (cond_split[0], "rep")
                else:
                    rep = None

                if "PreProcessReads" in job:
                    print("in preprocess reads")
                    path = get_exp_info(mirna, exp, cond, nb=nb, tp=tp,
                                        rep=rep)[0]
                else:
                    if "-uniq" in job:
                        read_dir = "reads_unique"
                    else:
                        read_dir = "reads"

                    path = get_analysis_path(mirna, exp, cond, analysis_type=read_dir)
                if "I_combined" in path:
                    path_update = path.split("I_combined.txt")[0] + "I.txt"
                elif "I_TGT_combined" in path:
                    path_update = path.split("I_TGT_combined.txt")[0] + "I.txt"
                elif "I_ACA_combined" in path:
                    path_update = path.split("I_ACA_combined.txt")[0] + "I.txt"
                elif "0_combined" in path:
                    path_update = path.split("0_combined.txt")[0] + "0.txt"
                else:
                    path_update = path
                total_lines = count_lines(path_update)
                host_chosen = False
                while not host_chosen:
                    print("getting joblist")
                    _joblist = JobList()
                    jobs_done = set(jobs_active.keys()).difference(set(_joblist.ids))
                    for job_id in jobs_done:
                        time_complete = time.time() - jobs_active[job_id][0]
                        threads = jobs_active[job_id][1]
                        job_lines = jobs_active[job_id][2]
                        del jobs_active[job_id]
                        jobs_completed.append((job_id, time_complete, job_lines, threads))
                    print("past job id loop")
                    print(len(jobs_done))
                    if len(jobs_done) > 0:
                        print("\n".join(["\t".join([str(j) for j in i])
                                         for i in jobs_completed]))
                    _hostlist = HostList()
                    print("Made host list object; about to assess which are active.")
                    len_active = len(_hostlist["active"])
                    print("number of active hosts: %s" %len_active)
                    if (len_active == 0):
                        print("sleeping")
                        time.sleep(0.1)
                    else:
                        print(311)
                        best_host = _hostlist["active"][0]
                        # if best_host == prior_host:
                        #     time.sleep(0.1)
                        #     print("Prior host:")
                        if int(best_host[1]) > 280 - len(_joblist):
                            print("max jobs")
                            time.sleep(0.1)
                        elif int(best_host[1]) == 1:
                            print("can't parallelize")
                            time.sleep(0.1)
                        else:
                            host_chosen = True

                print("host chosen")
                n_jobs = min(best_host[1], 21)
                bsub_line = "bsub -m %s -n %s %s -jobs %s" %(best_host[0], n_jobs, job,
                                                             n_jobs - 1)
                prior_host = best_host
                job_output = Popen(bsub_line, shell=True, stdout=PIPE).communicate()[0]
                if "is submitted to default queue <normal>" in job_output:
                    print("job submitted")
                    job_id = job_output.split("<")[1].split(">")[0]
                    jobs_active[job_id] = (time.time(), best_host[1], total_lines)
                    jobs.remove(job)
                    time.sleep(3)
            else:
                if say_finished:
                    print("DONE SUBMITTING JOBS")
                    say_finished = False
                _joblist = JobList()
                jobs_done = set(jobs_active.keys()).difference(set(_joblist.ids))
                for job_id in jobs_done:
                    time_complete = time.time() - jobs_active[job_id][0]
                    threads = jobs_active[job_id][1]
                    job_lines = jobs_active[job_id][2]
                    del jobs_active[job_id]
                    jobs_completed.append((job_id, time_complete, job_lines, threads))
                if len(jobs_done) > 0:
                    print("\n".join(["\t".join([str(j) for j in i])
                                     for i in jobs_completed]))
                    print("one minute")
                    time.sleep(60)


    os.remove(master_queue_path)
    # with open("temp_out_%s.txt" %(time.time()), "w+") as file_out:
    #     file_out.write("\n".join(["\t".join([str(j) for j in i])
    #                                  for i in jobs_completed]))
if __name__ == "__main__":
    main()

