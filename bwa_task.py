import os
import subprocess
from datetime import datetime
import genomon_rc as res


class Bwa_task():

    self.task_name = "bwa_mem"
    self.script_resource = """

#!/bin/bash
#
# Set SGE
#
#$ -S /bin/bash         # set shell in UGE
#$ -cwd                 # execute at the submitted dir

{bwa} mem -T {min_score} {ref_fa} {fastq1} {fastq2}  > {sam}

    """

class Stage_task():

    def __init__(self, qsub_option, script_dir):
       self.qsub_option = qsub_option
       self.script_dir = script_dir


    def ruffus_run(self, **kwargs):
        # Make shell script

        self.script_resource = self.get_script_resource(self):
        now = datetime.now()
        shell_script_name = res.file_timestamp_format.format(
                                 name=self.task_name,
                                 year=now.year,
                                 month=now.month,
                                 day=now.day,
                                 hour=now.hour,
                                 min=now.minute,
                                 msecond=now.microsecond )
        
        shell_script_full_path = "{script}/{file}.sh".format(script = self.script_dir, file = shell_script_name)
        shell_script_file = open( shell_script_full_path, 'w' )
        shell_script_file.write( self.script_resource.format( **kwargs))
        shell_script_file.close()

        shell_script = self.get_script_resource(self):
        process = subprocess.Popen( 'qsub  -sync yes -now no {qsub_option} ' + shell_script_file_path,
                                    qsub_option=self.qsub_option,
                                    shell=True,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE )

        std_out, std_err = process.communicate()

        print process.returncode 
    
