""" Tools to use htcondor via python.
Requires you to be on a computer where htcondor is set up.
For local machines see: https://twiki.cern.ch/twiki/bin/view/ABPComputing/LxbatchHTCondor

Hint: I had to add the path of the newly installed modules to the batch_Krb5_credential script to
be able to find the Authen::Krb5 package:
(use lib "/home/jdilly/perl5/lib/perl5/x86_64-linux-gnu-thread-multi/";)

Also, make sure the SCHEDD_NAME is set properly in the config file mentioned above.
(i.e. go to lxplus, run condor_q and see which scheduler was assigned to you)

Python Readme:
https://htcondor-python.readthedocs.io/en/latest

IMPORTANT: This functionality relies on shared-space of all files between htcondor and the user.
"""
import subprocess

import madx_wrapper
import os
import htcondor
import six
from utils import logging_tools
from utils.contexts import suppress_exception

LOG = logging_tools.get_logger(__name__)

SHEBANG = "#!/bin/bash"
BASHFILE_MASK = "htcbash.{:s}.sh"
SUBFILE = "queuehtc.sub"
MADX_PATH = madx_wrapper.MADX_AFS_PATH  # always use afs

CMD_SUBMIT = "condor_submit"
CMD_VIEW = "condor_q"
CMD_RM = "condor_rm"


# HTCondor Python Methods ######################################################


def submit_job(job):
    """ Submits the job to the scheduler
    This would be the nicer way, BUT: the python bindings do currently not forward the
    keberos token.
    """
    raise NotImplementedError("Submitting via python not implemented because of kerberos")
    schedd = htcondor.Schedd()  # works if schedd_name is set up
    with schedd.transaction() as txn:
        clusterId = job.queue(txn)  # assuming queue args are set in job, or 1
    return clusterId


def view_history():
    schedd = htcondor.Schedd()
    for ad in schedd.history('true', ['Owner', 'ProcId', 'ClusterId', 'JobStatus', 'WallDuration'], 2):
       LOG.info(ad)


# Subprocess Methods ###########################################################


def create_subfile_from_job(folder, job):
    """ Write file to submit to htcondor """
    subfile = os.path.join(folder, SUBFILE)

    with open(subfile, "w") as f:
        f.write(str(job))
    return subfile


def submit_jobfile(jobfile):
    """ Submit subfile to htcondor via subprocess """
    LOG.info("Sending {:s} to htcondor.".format(jobfile))
    _start_subprocess_with_logger([CMD_SUBMIT, jobfile])


def log_current_jobs_subprocess():
    _start_subprocess_with_logger(CMD_VIEW)


# Job Creation #################################################################


def create_multijob_for_bashfiles(folder, n_files, duration="longlunch"):
    """ Function to create a HTCondor job assuming n_files bash-files. """
    dura_key, dura_val = _get_duration(duration)

    job = htcondor.Submit({
        "MyId": "htcondor",
        "executable": os.path.join(folder, BASHFILE_MASK.format("$(ProcId)")),
        "arguments": "$(ClusterId) $(ProcId)",
        "initialdir": os.path.join(folder),
        "output": os.path.join("$(initialdir)", "$(MyId).$(ClusterId).$(ProcId).out"),
        "error": os.path.join("$(initialdir)", "$(MyId).$(ClusterId).$(ProcId).err"),
        "log": os.path.join("$(initialdir)", "$(MyId).$(ClusterId).$(ProcId).log"),
        dura_key: dura_val,
    })
    job.setQArgs("{:d}".format(n_files))
    return job


def create_job_for_bashfile(bashfile, duration="longlunch"):
    """ Returns a simple Submit() object for the bashfile. """
    dura_key, dura_val = _get_duration(duration)
    folder = os.path.dirname(bashfile)
    job = htcondor.Submit({
        "MyId": os.path.splitext(os.path.basename(bashfile))[0],
        "executable": bashfile,
        "arguments": "$(ClusterId)",
        "initialdir": os.path.join(folder),
        "output": os.path.join("$(initialdir)", "$(MyId).$(ClusterId).out"),
        "error": os.path.join("$(initialdir)", "$(MyId).$(ClusterId).err"),
        "log": os.path.join("$(initialdir)", "$(MyId).$(ClusterId).log"),
        dura_key: dura_val,
    })
    job.setQArgs("{:d}".format(1))
    return job


# For MADX #####################################################################


def write_madx_bash(folder, id, madx_files):
    """ Write bash file to call madx files

    Args:
        folder: Folder to write file into
        id: id for current bash_file
        madx_files: List of madx-files to run in this job.

    Returns:
        Path to current madx file

    """
    with suppress_exception(TypeError):
        id = "{:02d}".format(id)  # in case it's an int

    filename = os.path.join(folder, BASHFILE_MASK.format(id))
    with open(filename, "w") as f:
        f.write(SHEBANG + "\n")
        for madx_file in madx_files:
            f.write("{:s} < {:s}\n".format(MADX_PATH, madx_file))
    return filename


# Helper #######################################################################


def _get_duration(duration):
    if isinstance(duration, six.string_types):
        # espresso 20min, microcentury 1h, longlunch 2h, workday 8h,
        # tomorrow 1d, testmatch 3d, nextweek 1w
        return "+JobFlavour", '"{:s}"'.format(duration)
    return "+MaxRuntime", str(duration)  # runtime in seconds


def _start_subprocess_with_logger(command):
    process = subprocess.Popen(command, shell=False,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT,)
    for line in iter(process.stdout.readline, b''):
        line = line.strip()
        if line:
            LOG.info(line)
    status = process.wait()
    return status


# Script Mode ##################################################################


if __name__ == '__main__':
    raise EnvironmentError("{:s} is not supposed to run as main.".format(__file__))
    # afs_folder = "/afs/cern.ch/work/j/jdilly/public/"
    # job = create_job_for_bashfile(afs_folder + "job.sh")
    # sub = create_subfile_from_job(afs_folder, job)
    # submit_jobfile(sub)
    # submit_job(job)
    # view_history()