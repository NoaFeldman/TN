
# switch terminal to defaul (English) - Wb,Dec29,09
  export LANG=en_US.UTF-8 # de_DE.UTF-8
  export LC_ALL=$LANG
# export LANGUAGE=$LANG
  unset LANGUAGE

  if [ $(echo $PATH | grep -c $HOME/bin) -eq 0 ]; then
     export PATH="$PATH:.:$HOME/bin"
  fi


# module no longer `exported' by default // Wb,Sep08,16
  if [ `typeset -F module` ]; then
     export -f module
  fi

# -------------------------------------------------------------------- #
# ALIASES AND FUNCTIONS
# -------------------------------------------------------------------- #

  alias ls='ls --color'    # add colors for filetype recognition
  alias la='lt -a'         # based on bin/lt

# LINES would not be exported to shell scripts
  alias clr='i=0; while [ $i -le $LINES ]; do printf "\n"; i=$[i+1]; done; clear'

# default using -p : --preserve=mode,ownership,timestamps
  alias cp='cp -i --preserve=mode,timestamps'
  alias rm='rm -i'
  alias mv='mv -i'
  alias path='echo -e "\n  ${PATH//:/\\n  }"'
  alias ldpath='echo -e "\n  ${LD_LIBRARY_PATH//:/\\n  }"'

  export PS1='\n\s-\A \u@\h/\w\n$ '

# -------------------------------------------------------------------- #
# -------------------------------------------------------------------- #

