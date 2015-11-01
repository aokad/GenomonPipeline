#! /usr/bin/env python

import datetime
import subprocess

file_timestamp_format = "{name}_{year:0>4d}{month:0>2d}{day:0>2d}_{hour:0>2d}{min:0>2d}_{msecond:0>6d}"

class Stage_task(object):

    def __init__(self, qsub_option, script_dir):
        self.qsub_option = qsub_option
        self.script_dir = script_dir


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
        shell_script_file = open(shell_script_full_path, 'w')
        shell_script_file.write(self.script_template.format(**arguments))
        shell_script_file.close()

        qsub_commands = ['qsub', '-sync', 'yes']
        if max_task != 0:
            qsub_commands.extend(['-t', '1-'+str(max_task)+':1'])

        qsub_options = self.qsub_option.split(' ')
        returncode = subprocess.call(qsub_commands + qsub_options + [shell_script_full_path])

        return returncode
        # if returncode != 0:
        #     raise


