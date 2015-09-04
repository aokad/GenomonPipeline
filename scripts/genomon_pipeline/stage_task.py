#! /usr/bin/env python

import datetime

file_timestamp_format = "{name}_{year:0>4d}{month:0>2d}{day:0>2d}_{hour:0>2d}{min:0>2d}_{msecond:0>6d}"

class Stage_task(object):

    def __init__(self, qsub_option, script_dir):
        self.qsub_option = qsub_option
        self.script_dir = script_dir


    def task_exec(self, arguments):
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

        """
        process = subprocess.Popen('qsub -sync yes -now no {qsub_option} ' + shell_script_file_path,
                                    qsub_option=self.qsub_option,
                                    shell=True,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)

        std_out, std_err = process.communicate()
    
        print process.returncode 
        """

