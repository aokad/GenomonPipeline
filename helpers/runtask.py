#!/usr/bin/python
"""

Genome Analysis Pipeline
RunTask class


"""

import subprocess
from resource import genomon_rc as rc

#
# Run process
#
class RunTask:

    def __init__(   self,
                    run_mode = None,
                    ncpus = 0,
                    log = None,
                    resubmit = False,
                    qsub_cmd = 'qsub',
                    drmaa_native = None,
                    log_dir = None,
                    work_dir = None
                ):
        """
        Constructor

        """

        self.run_mode = run_mode
        self.log = log
        self.resubmit = resubmit
        self.qsub_cmd = qsub_cmd

        self.log.info( "RunTask: run_mode={run_mode}".format( run_mode = run_mode ) )

        if run_mode == 'MPI':
            self.ncpus = ncpus
            self.comm = []

            from mpi4py import MPI
            self.MPI = MPI

            self.log.info( "RunTask: ncpus={cpus}".format( ncpus = ncpus ) )

        elif run_mode == 'DRMAA':
            from drmaa_manage import JobManage as JM
            self.drmaa= JM( native_param = drmaa_native, log_dir = log_dir, work_dir = work_dir )

    def __del__( self ):
        """
        Destructor

        """
        if self.run_mode == 'MPI':
            self.disconnect()
        elif self.run_mode == 'DRMAA':
            self.drmaa.delete_job_template()


    def run_arrayjob( self, run_cmd, cmd_options, id_start = 1, id_end = 1, id_step = 1 ):
        return_code = 0
        if self.run_mode == 'MPI':
            return_code = self.__runtask_by_mpi(
                                run_cmd,
                                id_start = id_start,
                                id_end = id_end,
                                id_step = id_step )
        elif self.run_mode == 'DRMAA':
            return_code = self.__runtask_by_drmaa(
                                run_cmd,
                                cmd_options,
                                id_start = id_start,
                                id_end = id_end,
                                id_step = id_step )
        else:
            run_cmd_tmp = "-t {id_start}-{id_end}:{id_step}, {cmd}".format(
                                id_start = id_start,
                                id_end = id_end,
                                id_step = id_step,
                                cmd = run_cmd )
            return_code = self.__runtask_by_qsub( run_cmd_tmp, cmd_options )

        return return_code

    def runtask( self, run_cmd, cmd_options ):
        """
        Front end funtion to run task

        """

        if self.run_mode == 'MPI':
            return_code = self.__runtask_by_mpi( run_cmd, cmd_options )
        elif self.run_mode == 'DRMAA':
            return_code = self.__runtask_by_drmaa( run_cmd, cmd_options )
        else:
            return_code = self.__runtask_by_qsub( run_cmd, cmd_options )

        return return_code

    def disconnect( self ):
        id = 0
        for comm_tmp in self.comm:
            self.log.info( "Disconnect process {id}".format( id = id ) )
            comm_tmp.Disconnect()
            id += 1

    def __runtask_by_mpi( self,
                          run_cmd,
                          id_start = 1,
                          id_end = 1,
                          id_step = 1 ):
        """
        Submit a job by MPI.

        """
        self.log.info( 'DRMAA ' )
        self.log.info( "command         = {cmd}".format(cmd = run_cmd) )

        return_code = 0

        try:

            if id_start < id_end:
                self.log.info( "bulk job\n" )
                for id in range( id_start, id_end, id_step ):
                    run_cmd = "SGE_TASK_ID={0}\n".format( id ) + run_cmd
                    self.comm.append( self.MPI.COMM_SELF.Spawn( sys.executable, args=['mpi_worker.py', run_cmd ], maxprocs=1 ) )
            else:
                self.log.info( "normal job\n" )
                self.comm.append( self.MPI.COMM_SELF.Spawn( sys.executable, args=['mpi_worker.py', run_cmd ], maxprocs=1 ) )

        except OSError as e:
            self.log.error( "RunTask.runtaskby_qsub failed." )
            self.log.error( "OS error." )
            return_code = 1

        except IOError as (errno, strerror):
            self.log.error( "RunTask.runtaskby_qsub failed." )
            self.log.error( "IOError {0}{1}",format( errno, strerror ) )
            return_code = 1

        except Exception as e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            self.log.error( "RunTask.runtaskby_qsub failed." )
            self.log.error( "Unexpected error." )
            self.log.error("{0}: {1}:{2}".format( exc_type, fname, exc_tb.tb_lineno) )
            return_code = 1

        return return_code

    def __runtask_by_drmaa( self,
                            run_cmd,
                            cmd_options,
                            id_start = 0,
                            id_end = 0,
                            id_step = 0 ):
        """
        Submit a job by DRMAA.

        """
        self.log.info( 'DRMAA ' )
        self.log.info( "command         = {cmd}".format(cmd = run_cmd) )
        self.log.info( "command options = {cmd_options}".format(cmd_options = cmd_options) )

        if id_step != 0:
            self.log.info( "bulk job: {0}-{1}:{2}\n".format( id_start, id_end, id_step ) )
            self.drmaa.run_array_job( run_cmd,
                                      cmd_options = cmd_options,
                                      id_start = id_start,
                                      id_end = id_end,
                                      id_step = id_step )
        else:
            self.log.info( "normal job\n" )
            self.drmaa.run_job( run_cmd,
                                cmd_options = cmd_options )

        return_value = self.drmaa.wait_jobs()

        return return_value

        

    def __runtask_by_qsub( self, run_cmd, cmd_options ):
        """
        Submit a job by qsub.

        """
        self.log.info( 'qsub ' )
        self.log.info( "command     = {cmd}".format(cmd = run_cmd) )
        self.log.info( "cmd_options = {cmd_options}\n".format(cmd_options = cmd_options) )

        p_return_code = 0
        while True:

            cmd_tmp = self.qsub_cmd.format(
                                cmd_options = cmd_options,
                                cmd     = run_cmd )

            std_out = None
            std_err = None
            try:
                process = subprocess.Popen( cmd_tmp,
                                            shell=True,
                                            stdout=subprocess.PIPE,
                                            stderr=subprocess.PIPE )

                std_out, std_err = process.communicate()
                p_return_code = process.returncode

                for return_str in std_out.split( '\n' )[1:]:
                    #
                    # -sync  y  causes  qsub to wait for the job to complete before exiting.  If the job completes 
                    # successfully, qsub's exit code will be that of the completed job.  If the job fails to complete
                    # suc cessfully,  qsub  will  print  out a error message indicating why the job failed
                    # and will have an exit code of 1.
                    # If qsub is interrupted, e.g. with CTRL-C, before the job completes, the job will be canceled.
                    #
                    # Special case:
                    #   Job was qdel-ed by user.
                    #   The process of 'qsub -sync y' returns 0, which is a normal exit.
                    #   We need to detect if job is qdel-ed.
                    #
                    if 0 == return_str.find( 'Unable to run job' ):
                        raise ValueError( 'Submitted job was terminated.' )

            except OSError as e:
                self.log.error( "RunTask.runtaskby_qsub failed." )
                self.log.error( "OS error." )
                p_return_code = 1

            except IOError as (errno, strerror):
                self.log.error( "RunTask.runtaskby_qsub failed." )
                self.log.error( "IOError: {0}{1}".format( errno, strerror ) )
                p_return_code = 1

            except ValueError as e:
                self.log.error( e )
                p_return_code = 1

            except:
                self.log.error( "RunTask.runtaskby_qsub failed." )
                self.log.error( "Unexpected error." )
                p_return_code = 1

            self.log.info( "STDOUT: {stdout}".format( stdout = std_out ) )
            self.log.info( "STDERR: {stderr}\n".format( stderr = std_err ) )

            if ( p_return_code == 0 or not self.resubmit ):
                break

        return p_return_code

