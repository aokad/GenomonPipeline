#! /usr/bin/env python

import os, datetime, subprocess
from genomon_pipeline.config.run_conf import *
import drmaa
 
file_timestamp_format = "{name}_{year:0>4d}{month:0>2d}{day:0>2d}_{hour:0>2d}{min:0>2d}_{msecond:0>6d}"

class Stage_task(object):

    def __init__(self, qsub_option, script_dir):
        self.qsub_option = qsub_option
        self.script_dir = script_dir
        self.drmaa = run_conf.drmaa


    def task_exec(self, arguments, max_task=0):
        # Make shell script

        now = datetime.datetime.now()
        shell_script_name = file_timestamp_format.format(
                                 name=self.task_name,
                                 year=now.year,
                                 month=now.month,
                                 day=now.day,
                                 hour=now.hour,
                                 min=now.minute,
                                 msecond=now.microsecond )

        shell_script_full_path = "{script}/{file}.sh".format(script = self.script_dir, file = shell_script_name)
        print shell_script_full_path
        shell_script_file = open(shell_script_full_path, 'w')
        shell_script_file.write(self.script_template.format(**arguments))
        shell_script_file.close()
        os.chmod(shell_script_full_path, 0750)

        if self.drmaa:

            s = drmaa.Session()
            s.initialize()

            jt = s.createJobTemplate()
            jt.jobName = shell_script_name
            jt.outputPath = ':' + run_conf.project_root + '/log'
            jt.errorPath = ':' + run_conf.project_root + '/log' 
            jt.nativeSpecification = self.qsub_option
            jt.remoteCommand = shell_script_full_path

            print jt.jobName
            print jt.outputPath
            print jt.nativeSpecification
            print jt.remoteCommand

            jobid = s.runJob(jt)
            now = datetime.datetime.now()
            date = now.strftime("%Y-%m-%d %H:%M")
            print "Date/Time: " + date 
            print "Job has been submitted with id: " + jobid + " at Date/Time: " + date
            retval = s.wait(jobid, drmaa.Session.TIMEOUT_WAIT_FOREVER)
            now = datetime.datetime.now()
            date = now.strftime("%Y-%m-%d %H:%M")
            print "Job: " + str(retval.jobId) + ' finished with status: ' + str(retval.hasExited) + ' and exit status: ' + str(retval.exitStatus) + " at Date/Time: " + date
            s.deleteJobTemplate(jt)
            s.exit()
            returncode = retval.exitStatus

        else:
            """
            qsub_commands = ['qsub', '-sync', 'yes']
            if max_task != 0:
                qsub_commands.extend(['-t', '1-'+str(max_task)+':1'])

            qsub_options = self.qsub_option.split(' ')
            returncode = subprocess.call(qsub_commands + qsub_options + [shell_script_full_path])
            """
            returncode = 1
        if returncode != 0:
            raise


