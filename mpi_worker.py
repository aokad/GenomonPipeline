import sys
from mpi4py import MPI
import subprocess

return_code = 0

try:
    print "mpi_worker: {0} started".format( sys.argv[ 1 ] )
    comm = MPI.Comm.Get_parent()
    p = subprocess.Popen( sys.argv[ 1 ],
                          shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          close_fds=True )
    stdout, stderr = p.communicate()
    print "mpi_worker:STDOUT: {0}".format( stdout )
    print "mpi_worker:STDERR: {0}".format( stderr )
    print "\n"
    comm.Disconnect()

except OSError as e:
    log.error( "mpi_worker failed." )
    log.error( "OS error." )
    return_code = 1

except IOError as (errno, strerror):
    log.error( "mpi_worker failed." )
    log.error( "IOError {0}{1}",format( errno, strerror ) )
    return_code = 1

except Exception as e:
    exc_type, exc_obj, exc_tb = sys.exc_info()
    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    log.error( "mpi_worker failed." )
    log.error( "Unexpected error." )
    log.error("{0}: {1}:{2}".format( exc_type, fname, exc_tb.tb_lineno) )
    return_code = 1

sys.exit( return_code )

