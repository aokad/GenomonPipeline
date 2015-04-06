#!/usr/bin/python
"""

Genome Analysis Pipeline
RunTask class


"""

import subprocess
import genomon_rc as rc

#
# Run process
#
class RunTask:

    def __init__(   self,
                    enable_mpi = False,
                    ncpus = 0,
                    log = None,
                    resubmit = False,
                    max_mem = 32,
                    qsub_cmd = None,
                ):
        """
        Constructor

        """
        self.enable_mpi = enable_mpi
        self.log = log
        self.resubmit = resubmit
        self.max_mem = max_mem
        self.qsub_cmd = qsub_cmd

        if enable_mpi:
            self.ncpus = ncpus
            self.comm = []

            from mpi4py import MPI
            self.MPI = MPI

            self.log.info( "RunTask: mpi={enable_mpi}".format( enable_mpi = enable_mpi ) )
            self.log.info( "RunTask: ncpus={cpus}".format( ncpus = ncpus ) )


    def __del__( self ):
        """
        Destructor

        """
        if self.enable_mpi:
            self.disconnect()


    def run_arrayjob( self, job_queue, memory, run_cmd, id_list ):
        return_code = 0
        if self.enable_mpi:
            pass
        else:
            run_cmd_tmp = "-t {id_list}, {cmd}".format(
                                id_list = id_list,
                                cmd = run_cmd )
            return_code = self.__runtask_by_qsub( job_queue, memory, run_cmd_tmp )

        return return_code

    def runtask( self, job_queue, memory, run_cmd ):
        """
        Front end funtion to run task

        """

        if self.enable_mpi:
            return_code = self.__runtask_by_mpi( job_queue, memory, run_cmd )
        else:
            return_code = self.__runtask_by_qsub( job_queue, memory, run_cmd )

        return return_code

    def disconnect( self ):
        id = 0
        for comm_tmp in self.comm:
            self.log.info( "Disconnect process {id}".format( id = id ) )
            comm_tmp.Disconnect()
            id += 1

    def __runtask_by_mpi( self, job_queue, memory, run_cmd ):

        return_code = 0

        try:
            self.comm.append( self.MPI.COMM_SELF.Spawn( sys.executable, args=['mpi_worker.py', run_cmd ], maxprocs=1 ) )

        except OSError as e:
            self.log.error( "RunTask.runtaskby_qsub failed." )
            self.log.error( "OS error." )
            return_code = 1

        except IOErr as (errno, strerror):
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

    def __runtask_by_qsub( self, job_queue, memory, run_cmd ):
        """
        Submit a job by qsub.

        """
        self.log.info( '# runtask_by_qsub' )
        self.log.info( "command = {cmd}".format(cmd = run_cmd) )
        self.log.info( "memory  = {mem}".format(mem = memory) )
        self.log.info( "job     = {job}\n".format(job = job_queue) )

        if job_queue == 'mjob':
            job_queue = ''

        return_code = 0
        p_return_code = 0
        while True:

            cmd_tmp = self.qsub_cmd.format(
                                s_vmem  = memory,
                                mem_req = memory[:-1],
                                job_queue = job_queue,
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

            except OSError as e:
                self.log.error( "RunTask.runtaskby_qsub failed." )
                self.log.error( "OS error." )
                return_code = 1

            except IOErr as (errno, strerror):
                self.log.error( "RunTask.runtaskby_qsub failed." )
                self.log.error( "IOError: {0}{1}".format( errno, strerror ) )
                return_code = 1

            except CalledProcessError as e:
                self.log.error( "CalledProcessError: return code: {id}".format( id = e.returncode ) )
                return_code = 1

            except:
                self.log.error( "RunTask.runtaskby_qsub failed." )
                self.log.error( "Unexpected error." )
                return_code = 1

            else:
                return_code = 0
                    
            self.log.info( "STDOUT: {stdout}".format( stdout = std_out ) )
            self.log.info( "STDERR: {stderr}".format( stderr = std_err ) )

            memory = str( int( memory[0:-1] ) * 2 ) + 'G'

            if ( return_code == 0 or
                 not self.resubmit or
                 int( memory[0:-1] ) > self.max_mem
               ):
                break

        return return_code

