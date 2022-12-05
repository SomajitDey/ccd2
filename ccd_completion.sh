# Source me for Tab based autocomplete on `ccd <Tab>`

#TODO: Extract the options list automatically from the codebase, i.e. from ccd_<option>.f90
_ccd_completions()
{
  COMPREPLY=($(compgen -W "run rinit traj_to_legacy cpt_to_xy status visual" "${COMP_WORDS[1]}"))
}

complete -F _ccd_completions ccd
