# Exit script for ccd
# This script is sourced by the ccd driver
# Best understood in conjuction with ccd_config.sh

# Remove the hidden parameter file actually used by the sim engine
rm -f .params.in

# Show quote
ccd_quote
