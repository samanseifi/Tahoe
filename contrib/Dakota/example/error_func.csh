#! /bin/csh

# $argv[1] is .in from Dakota
# $argv[2] is .out returned to Dakota

#echo $PATH

#---------------
# configuration 
#---------------
set TAHOE          = /Users/vicky/Tahoe/tahoe-install/tahoe/tahoe

set DAKOTA_UTIL    = /Users/vicky/Tahoe/tahoe-install/contrib/Dakota
set EXTRACT        = extract_1D
set TAHOE_TEMPLATE = ${DAKOTA_UTIL}/example/mat.2.a.xml.tmpl
set RESULT_EXT     = io1.run
set ERROR_FUNC     = least_square
set REF_RESULTS    = ${DAKOTA_UTIL}/example/fvsd

#--------------
# create tahoe input file
#--------------

perl ${DAKOTA_UTIL}/dakota2tahoe.pl $argv[1] ${TAHOE_TEMPLATE} $argv[1].xml

#-----
# run 
#-----

${TAHOE} -f $argv[1].xml > $argv[1].console

#---------------
# extract data
#---------------

${DAKOTA_UTIL}/${EXTRACT} $argv[1].${RESULT_EXT} $argv[1].${RESULT_EXT} $argv[1]_Z.fvsd

#---------------
# compute error function
#---------------

${DAKOTA_UTIL}/${ERROR_FUNC} ${REF_RESULTS} $argv[1]_Z.fvsd $argv[2]

#---------------------------------------------------------------
# clean up
#---------------------------------------------------------------

rm -f $argv[1].*xml $argv[1].io* $argv[1].out $argv[1].console $argv[1]*.fvsd  $argv[1]*.run.* >& /dev/null
