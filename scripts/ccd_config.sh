# Configuration script for ccd
# This script is sourced by the ccd driver
# Also see ccd_exit.sh

[[ -n "${CCD_PARAMS_PATH}" ]] || CCD_PARAMS_PATH="params.in"
[[ -n "${CCD_RC_PATH}" ]] || CCD_RC_PATH="${HOME}/.ccdrc"

export CCD_PARAMS_PATH CCD_RC_PATH

# Prepare the parameter file to be actually used by the sim engine
sed '1i &params ! Namelist header' <(cat "${CCD_RC_PATH}" "${CCD_PARAMS_PATH}" 2>/dev/null; echo) > .params.in
echo '/ ! Namelist footer' >> .params.in
