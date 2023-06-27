# Simply source me with `source setup_test_env.sh` or `. setup_test_env.sh`

whereami(){
  # Returns the absolute path of the directory this script is in
  case "${BASH_SOURCE}" in
    /*)
      echo -n "${BASH_SOURCE%/*}" ;;
    */*)
      echo -n "${PWD}/${BASH_SOURCE%/*}" ;;
    *)
      echo -n "${PWD}";;
  esac
}
this_script_is_at="$(whereami)"

export PATH="${this_script_is_at}/build:${this_script_is_at}/scripts:${this_script_is_at}/hooks:${this_script_is_at}/tests:${PATH}"

. "${this_script_is_at}/ccd_completion.sh"
ulimit -s unlimited
export OMP_STACKSIZE=500m

export CCD_NO_QUOTES=set

echo "Done. To avail --help|-h messages, install properly with: make && make install"
