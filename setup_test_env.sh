# Simply source me with `source setup_test_env.sh` or `. setup_test_env.sh`

export PATH="${PWD}/build:${PWD}/scripts:${PWD}/hooks:${PATH}"
. ccd_completion.sh
ulimit -s unlimited
echo "Done"
