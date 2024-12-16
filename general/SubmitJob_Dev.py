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
          stdout=PIPE).communicate()[0].split("\n"))[1:-1]

def count_lines(path):
    return(int(Popen("wc -l %s" %(path), shell=True,
          stdout=PIPE).communicate()[0].split()[0]))

def ls(path):
    return(Popen("ls %s" %(path), shell=True,
          stdout=PIPE, stderr=PIPE).communicate()[0].strip().split())

def check_jobs():
    return(ls("job_queue/*_jobs.txt"))

def check_master():
    return(ls("job_queue/*_master.txt"))

def get_job_ids():
    return(Popen("bjobs | cut -d ' ' -f 1", shell=True,
          stdout=PIPE).communicate()[0].split("\n"))[1:-1]

def get_job_hosts():
    return(Popen("bjobs | cut -c 46- | cut -d ' ' -f 1", shell=True,
          stdout=PIPE).communicate()[0].split("\n"))[1:-1]


def get_job_statuses():
    return(Popen("bjobs | cut -c 17- | cut -d ' ' -f 1", shell=True,
          stdout=PIPE).communicate()[0].split("\n"))[1:-1]


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
            print("here")
            return(self.container[key])
        else:
            return
    def __len__(self):
        return(sum([len(_job) for _job in self.jobs]))



def get_all_hosts():
    return([i.split()[0] for i in Popen("bhosts", shell=True,
                 stdout=PIPE).communicate()[0].strip().split("\n")[1:]])

def get_job_info(name):
    titles = Popen("bhosts | head -1", shell=True,
                 stdout=PIPE).communicate()[0].strip().split()
    host_info = Popen("bhosts | grep %s" %(name), shell=True,
                 stdout=PIPE).communicate()[0].strip().split()
    return {key : value for key, value in zip(titles, host_info)[1:]}



class Host():
    def __init__(self, name):
        self.name = name
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
    EXCLUDED_HOSTS = ["slx-c03b04", "slx-c03b05", "slx-c03b06"]
    def __init__(self):
        self.list = get_all_hosts()
        self.hosts = {name : Host(name) for name in get_all_hosts()}
    def __len__(self):
        return(len([name for name in self.list
                    if self.hosts[name]["status"] == "ok"
                    and self.hosts[name]["jobs"] >= 0
                    and name not in HostList.EXCLUDED_HOSTS]))
    def __getitem__(self, key):
        if key == "active":
            return(sorted([(name, self.hosts[name]["jobs"], self.hosts[name]["status"])
                           for name in self.list
                           if self.hosts[name]["status"] == "ok"
                           and self.hosts[name]["jobs"] >= 0
                           and name not in HostList.EXCLUDED_HOSTS],
                          key = lambda host: host[1], reverse=True))

def main():
    time_start = time.time()
    arguments = ["job", "-test_job_binary"]
    args = parse_arguments(arguments)
    jobs = args[0].strip().split("\\n")
    test_job = args[1]

    if jobs[-1] == "":
        jobs = jobs[:-1]
    script_id = "%s_%s" %(time.time(), random.random())
    master_queue_path = "job_queue/%s_master.txt" %(script_id)
    job_queue_path = "job_queue/%s_jobs.txt" %(script_id)
    if check_master() or test_job:
        print("giving jobs")
        with open(job_queue_path, "w+") as queue_file:
            queue_file.write("\n".join(jobs))
        print(job_queue_path)
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
        tick = 0
        prior_host = ""
        host = ""
        jobs_accumulated = 0
        jobs_completed = []
        jobs_active = {}
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
                time.sleep(2)
                outside_jobs = check_jobs()
                print(outside_jobs)

            if len(jobs) > 0:
                with open(master_queue_path, "w+") as master_file:
                    master_file.write("\n".join(jobs))
                # print(len(jobs))
                job = jobs[0]
                mirna, exp, cond = job.split()[2:5]
                path = get_analysis_path(mirna, exp, cond, analysis_type="reads")
                total_lines = count_lines(path)
                host_chosen = False
                while not host_chosen:
                    _joblist = JobList()
                    jobs_done = set(jobs_active.keys()).difference(set(_joblist.ids))
                    for job_id in jobs_done:
                        # print("finished:")
                        # print(jobs_active[job_id])
                        # print("new time:")
                        # print(time.time())
                        time_complete = time.time() - jobs_active[job_id][0]
                        threads = jobs_active[job_id][1]
                        job_lines = jobs_active[job_id][2]
                        del jobs_active[job_id]
                        jobs_completed.append((job_id, time_complete, job_lines, threads))
                    if len(jobs_done) > 0:
                        print("\n".join(["\t".join([str(j) for j in i])
                                         for i in jobs_completed]))
                    _hostlist = HostList()
                    # print(len(_hostlist["active"]))
                    if (len(_hostlist["active"]) == 0):
                        # print("sleep no hosts")
                        time.sleep(0.1)
                    else:
                        best_host = _hostlist["active"][0]
                        if best_host == prior_host:
                            time.sleep(0.1)
                            print("Prior host:")
                        elif int(best_host[1]) > 280 - len(_joblist):
                            print("max jobs")
                            time.sleep(0.1)
                        elif int(best_host[1]) == 1:
                            print("can't parallelize")
                            time.sleep(0.1)
                        else:
                            host_chosen = True

                print(best_host)
                n_jobs = best_host[1]
                # n_jobs = max(20, n_jobs)
                bsub_line = "bsub -m %s -n %s %s -jobs %s" %(best_host[0], n_jobs, job,
                                                             n_jobs - 1)
                prior_host = best_host
                # print(bsub_line)
                job_output = Popen(bsub_line, shell=True, stdout=PIPE).communicate()[0]
                # print(job_output)
                if "is submitted to default queue <normal>" in job_output:
                    job_id = job_output.split("<")[1].split(">")[0]
                    # print(job_id)
                    jobs_active[job_id] = (time.time(), best_host[1], total_lines)
                    print("time_start:")
                    print(jobs_active[job_id][0])
                    print("total jobs:")
                    print(n_jobs)
                    print("total_lines:")
                    print(total_lines)
                    jobs.remove(job)
                    time.sleep(3)
            else:
                if say_finished:
                    print("DONE SUBMITTING JOBS")
                    say_finished = False


                _joblist = JobList()
                jobs_done = set(jobs_active.keys()).difference(set(_joblist.ids))
                for job_id in jobs_done:
                    print("finished:")
                    print(jobs_active[job_id])
                    print("new time:")
                    print(time.time())
                    time_complete = time.time() - jobs_active[job_id][0]
                    threads = jobs_active[job_id][1]
                    job_lines = jobs_active[job_id][2]
                    print("time_done:")
                    print(time_complete)
                    print("time taken:")
                    print(time_complete - jobs_active[job_id][0])
                    print("total lines:")
                    print(job_lines)
                    print("job threads:")
                    print(threads)
                    del jobs_active[job_id]
                    jobs_completed.append((job_id, time_complete, job_lines, threads))
                if len(jobs_done) > 0:
                    print("\n".join(["\t".join([str(j) for j in i])
                                     for i in jobs_completed]))
                    time.sleep(1)

            tick += 1
        os.remove(master_queue_path)
        with open("temp_out_%s.txt" %(time.time()), "w+") as file_out:
            file_out.write("\n".join(["\t".join([str(j) for j in i])
                                         for i in jobs_completed]))
if __name__ == "__main__":
    main()

