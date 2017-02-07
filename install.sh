#!/bin/bash

#   Mega2: Manipulation Environment for Genetic Analysis
#   Copyright (C) 1999-2016 Robert Baron, Justin R. Stickel, Charles P. Kollar,
#   Nandita Mukhopadhyay, Lee Almasy, Mark Schroeder, William P. Mulvihill,
#   Daniel E. Weeks, and University of Pittsburgh
#  
#   This file is part of the Mega2 program, which is free software; you
#   can redistribute it and/or modify it under the terms of the GNU
#   General Public License as published by the Free Software Foundation;
#   either version 3 of the License, or (at your option) any later
#   version.
#  
#   Mega2 is distributed in the hope that it will be useful, but WITHOUT
#   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#   for more details.
#  
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#  
#   For further information contact:
#       Daniel E. Weeks
#       e-mail: weeks@pitt.edu
# 
# ===========================================================================

# NOTE: Please refer to 'srcdir/version.h' which also encodes the Mega2 version
# information.
#VERSION=v4.9.1

VERSION_H='./srcdir/version.h'
MEGA2VER=`grep MEGA2VER ${VERSION_H} | sed -n 's/^#define MEGA2VER \(.*\)/\1/p'`
MEGA2REV=`grep MEGA2REV ${VERSION_H} | sed -n 's/^#define MEGA2REV \(.*\)/\1/p'`
MEGA2PATCH=`grep MEGA2PATCH ${VERSION_H} | sed -n 's/^#define MEGA2PATCH \(.*\)/\1/p'`
# Get the current Mega2 version numbers from the source...
VERSION="v${MEGA2VER}.${MEGA2REV}.${MEGA2PATCH}"

# Minimum acceptable version of nplplot...
MINIMUM_REQUIRED_NPLPLOT_VERSION='4.5'

#
# shell script to install mega2

#
# The user can bypass certain questions asked by the script by setting the
# appropriate environment variables. So, if the user always wants the behavior
# associated with these variables they can be set one time in the user's shell
# RC file. The variables are:
#
# MEGA2_BINARY_DIRECTORY (if defined) should be set to the path that the Mega2
# binaries will be installed. If this is not a path (directory) an error will
# be given, and the script will exit.
#
# MEGA2_OVERWRITE_BINARIES (if defined) should be set to any character string,
# it will cause the install script to always overwrite the binaries.

function verify() {
    echo

    echo Checking for python
    a=( `$which python 2>/dev/null` )
    if [[ $a == no || ! $a ]]; then 
        echo WARNING: Could not find python
#        read -p "Enter full path to python program: " -a python_path 
        echo -n "Enter full path to python program: " 
        read -a python_path 
    else
        python_path=$a 
        echo OK: Found $a
    fi

    echo Checking for gawk/nawk/awk
    a=( `$which gawk 2>/dev/null` )
    if [[ $a  == no || ! $a ]]; then
    #    echo gawk not found, looking for nawk
        a=( `$which nawk 2>/dev/null` )
        if [[ $a == no || ! $a ]]; then
    #	echo gawk, nawk not found, looking for awk
            a=( `$which awk 2>/dev/null` )
            if [[ $a == no || ! $a ]]; then
                echo WARNING: Could not find gawk, awk or nawk ..
                echo "         Some Mega2-generated scripts may fail to run."
            else
                echo OK: Found $a
            fi
        else
            echo OK: Found $a
        fi
    else
        echo OK: Found $a
    fi

    echo Checking for perl
    a=( `$which perl 2>/dev/null` )
    if [[ $a == no || ! $a ]]; then
        echo WARNING: Could not find perl
#       read -p "Enter full path to perl executable: " -a perl_path
        echo -n "Enter full path to perl executable: "
        read -a perl_path
    else
        perl_path=$a
        echo OK: Found $perl_path
    fi

    echo Checking for R
    a=( `$which R 2>/dev/null` )
    if [[ $a == no || ! $a ]]; then
       echo WARNING: Could not find R
    elif [[ ${SAVE:-""} == "" ]]; then
        echo OK: Found $a
        echo Checking for nplplot library installed in R
        rm -f nplplot_installed_version.txt
        R --no-save 1>/dev/null 2>&1 <<RSCRIPT
sink("nplplot_installed_version.txt")
packageVersion("nplplot")
RSCRIPT
        if [[ -s "nplplot_installed_version.txt" ]] ; then
            NPLPLOT_INSTALLED_VERSION=`sed -n 's/^\[1\] \(.*\)/\1/p' nplplot_installed_version.txt`
            # The quotes that R produces are the all too cute ASCII \xE2, rather than something like '`
            NPLPLOT_INSTALLED_VERSION=`echo $NPLPLOT_INSTALLED_VERSION |  sed -n 's/.$//p' | sed -n 's/^.//p'`
            if [[ $NPLPLOT_INSTALLED_VERSION < $MINIMUM_REQUIRED_NPLPLOT_VERSION ]] ; then
                echo "ERROR: Found nplplot version $NPLPLOT_INSTALLED_VERSION"
                echo "ERROR: At least version $MINIMUM_REQUIRED_NPLPLOT_VERSION is required."
                echo "ERROR: Please start R and install nplplot using the following command."
                echo "ERROR: install.packages('nplplot')"
                exit
            fi
            echo "OK: Found nplplot version $NPLPLOT_INSTALLED_VERSION"
        else
            echo "WARNING: nplplot not installed. To install it you"
            echo "WARNING: will need to start R and run the command:"
            echo "WARNING: install.packages('nplplot')"
            echo "ERROR: Mega2 installation not completed - retry after installing 'nplplot'"
            exit
        fi
    fi

    echo Checking for csh
    a=`$which csh 2>/dev/null`
    if [[ $a == no || ! $a ]]; then
       echo WARNING: Could not find csh
    else
       echo OK: Found $a
    fi

    # if this environment variable is set to a directory, then just install there
    # and don't query the user
    if [[ $MEGA2_BINARY_DIRECTORY ]] ; then
        mega2_path=$MEGA2_BINARY_DIRECTORY
        if [ ! -d $mega2_path ]; then
            echo "The directory '$mega2_path'"
            echo "associated with the environment variable MEGA2_BINARY_DIRECTORY does not exist."
            echo "Please create the directory or unset the environment variable and run install.sh"
            echo "(this script) again."
            exit
        fi
        echo "Using path $mega2_path from the environment variable MEGA2_BINARY_DIRECTORY"
        echo "for $module to be installed."
    else 
        echo
        echo "Enter full path for $module to be installed."
        echo "You will need write access."
        echo -n "[default /usr/local/bin]: "
        read -a mega2_path
    fi

    if [[ $mega2_path == "" ]]; then
        mega2_path=/usr/local/bin
    fi

    if [ ! -d $mega2_path ]; then
        echo "The path '$mega2_path' is not a directory.  Please create the directory and"
        echo "and run install.sh (this script) again."
        exit
    fi

    if [[ "`echo $PATH | grep $mega2_path`" == "" ]] ; then
        echo "*** The directory $mega2_path is not in your PATH, you will need to add it."
        echo
    fi
}

if [[ $1 ]]; then
    scripts=$1;
else 
    scripts=all;
fi

if [[ $scripts == "all" ]]; then
    module="mega2 and helper scripts"
    scripts="mega2 scripts";

elif [[ $scripts == "mega2" || $scripts = "mega2compile" ]]; then
    module=mega2

elif  [[ $scripts == "scripts" ]]; then
    module="perl scripts"

else

    echo ERROR: Unrecognized component $scripts to install.sh ;
    echo "       run again in one of these three ways:";
    echo "       ./install.sh all";
    echo "       ./install.sh mega2";
    echo "       ./install.sh mega2compile";
    echo "       ./install.sh scripts";
    exit -1;
fi


# Look for perl in /usr/bin or /usr/local/bin first
# 

perl_scripts="Rsimwalk2.pl Rmerlin_regress.pl Rallegro.pl \
	    make_gen_table.pl make_hwe_table.pl \
	    merlin2sw2.pl mega2log2html.pl"; 

python_scripts="l2a.py Rmerlin.py";

OSTYPE=`uname -s`
if [[ $OSTYPE == Darwin ]]; then 
    OSTYPE=darwin
elif [[ $OSTYPE == Linux ]]; then
    OSTYPE=linux
elif [[ ${OSTYPE:0:6} == CYGWIN ]]; then
    OSTYPE=cygwin
elif [[ ${OSTYPE:0:10} == MINGW64_NT ]]; then
    OSTYPE=msys2
elif [[ ${OSTYPE:0:10} == MINGW32_NT ]]; then
    OSTYPE=mingw
elif [[ ${OSTYPE:0:7} == MSYS_NT ]]; then
    OSTYPE=msys2
fi
export OSTYPE
echo OS $OSTYPE
echo

which=which
pgm=none

# Try to figure out which flavor of Linux this is...
function pgm_linux() {
    v=`sed -n /etc/issue -e "1s/Ubuntu \\([0-9\\.]*\\) .*/\\1/p"`
    if [ "$v" != "" ] ; then
        pgm=mega2_${VERSION}_linux_ubuntu.${v}x${x}
        return
    fi
    if [ -f /etc/centos-release ] ; then 
        v=`sed -n /etc/centos-release -e "1s/^CentOS Linux release \\([0-9\\.]*\\) (Core)/\\1/p"`
        if [ "$v" != "" ] ; then
            pgm=mega2_${VERSION}_linux_centos.${v}
            return
        fi
    fi
    pgm=mega2_${VERSION}_linux
}

# export so that makefiles in subdirectories have access to these...
export MAKE=make
export CC=cc
for scr in $scripts; do
    if [[ $scr == "mega2" || $scr == "mega2compile" ]]; then
        if [[ $OSTYPE == darwin ]]; then 
#            v=`ls /Developer/SDKs/|sed -n -e "\\$s/MacOSX\\(.*\\).sdk/\\1/p"`
            # /Developer no longer exist in Xcode 4.3 so we need to get the OS number a different way...
	    v=`sw_vers | grep '^ProductVersion:' | awk 'BEGIN{FS=" "}{print $2}' | awk 'BEGIN{FS="."}{print $1 "." $2}'`
            pgm=mega2_${VERSION}_darwin.${v}
        elif [[ $OSTYPE == linux ]]; then
            if [ -f /etc/issue ]; then
                x=`file /sbin/init| sed -n - -e "1s/.*ELF \\([0-9]*\\)-bit.*/\\1/p"`
                pgm_linux
            fi
        elif [[ $OSTYPE == SunOS ]]; then
            # uname -v prints 'Generic_147148-26' on Solaris 10.1, and '11.1' on Solaris 11.1
            # uname -r prints '5.10' on Solaris 10.1, and '5.11' on Solaris 11.1
            v=`uname -v`
            x=`uname -i`
            pgm=mega2_${VERSION}_solaris.${v}_${x} # e.g. mega2_v4.5.9_solaris.11.1_i86pc
            export MAKE=gmake
            export CC=gcc
        elif [[ $OSTYPE == cygwin ]]; then
            v=`uname -r|sed -n -e "s/\\(.*\\)(.*)/\\1/p"`
            pgm=mega2_${VERSION}_cygwin.${v}
        elif [[ $OSTYPE == mingw ]]; then
            which=where
            v=`uname -r|sed -n -e "s/\\(.*\\)(.*)/\\1/p"`
            pgm=mega2_${VERSION}_mingw.${v}
        elif [[ $OSTYPE == msys2 ]]; then
            CC=gcc
            which=where
            v=`uname -r|sed -n -e "s/\\(.*\\)(.*)/\\1/p"`
            pgm=mega2_${VERSION}_msys2.${v}
        fi

        if [[ -x mega2_bin/$pgm && $pgm != "none" && $scr != "mega2compile" ]]; then
            mega2prog=mega2_bin/$pgm
            echo Not recompiling mega2.  Using $mega2prog

        elif [ -f srcdir/Makedefs_${OSTYPE}.template ]; then
            cd srcdir
#           cp "Makedefs_"${OSTYPE}".template" Makedefs_${OSTYPE}
#            echo compiling mega2
            echo "Attempting to compile mega2 because no pre-compiled binary was found for your operating system."
	    compile_prerequisites=present

	    a=( `$which $MAKE 2>/dev/null` )
            # Check to see if we have any make...
            if [[ $a == no || ! $a ]]; then 
                compile_prerequisites=missing
            else
                a=( `$MAKE --version 2>/dev/null | head -1` )
                # Check to see if we have GNU make...
		if [[ $a != "GNU" ]]; then
                    compile_prerequisites=missing
                fi
            fi
            if [[ $compile_prerequisites != "present" ]]; then
                echo
                echo "GNU make is required for compiling mega2."
                if [[ $OSTYPE == linux ]]; then
                    echo "On Red Hat/CentOS/Fedora/Oracle 5 systems:"
                    echo "sudo yum install make"
                    echo "on Debian systems:"
                    echo "apt-get install build-essential"
                elif [[ $OSTYPE == darwin ]]; then 
                    echo "The Xcode command line tools must be installed to compile mega2."
                else
                    echo "For more information on GNU make please see 'http://www.gnu.org/software/make/'"
                fi
	    fi

	    a=( `$which $CC 2>/dev/null` )
            # Check to see if there is a c compiler...
            if [[ $a == no || ! $a ]]; then 
                echo
                echo "A c compiler is required for compiling mega2."
                compile_prerequisites="missing"
                if [[ $OSTYPE == linux ]]; then
                    echo "On Red Hat/CentOS/Fedora/Oracle 5 systems:"
                    echo "sudo yum install gcc"
                    echo "on Debian systems:"
                    echo "apt-get install build-essential"
                elif [[ $OSTYPE == darwin ]]; then 
		    a=( `$which xcode-select 2>/dev/null` )
		    # Check to see if Xcode is installed...
		    if [[ $a == no || ! $a ]]; then 
			echo "It appears that you do not have Xcode installed to allow you to compile mega2."
			echo "For the latest version of Xcode the command line tools are optional."
			echo "Install them using the Components tab of the Downloads preference panel of Xcode."
#		    elif [[ -n "`xcode-select -print-path | grep /Applications/Xcode.app`" ]]; then
		    elif [[ -d /Applications/Xcode.app ]]; then
			echo "For the current version of Xcode the 'Command Line Tools' are optional."
			echo "Install them from 'https://developer.apple.com/downloads/index.action'"
		    else
			echo "It appears that a version of Xcode is installed but the c compiler is not."
			echo "The command line tools must be installed to compile mega2."
		    fi
	        fi
	    fi

	    if [[ $compile_prerequisites != "present" ]]; then
		echo
		echo "See the Mega2 documentation for details on how to compile from source."
		exit
            fi

            if [[ ${SAVE:-""} == "" ]]; then
                $MAKE all
            elif [[ ${SAVE:-""} == "static" ]]; then
                $MAKE LN="" STATIC="use" distclean all
            else
                $MAKE LN="" distclean all
            fi
            make_status=$?
            if [[ $make_status > 0 ]]; then
                echo make program exit value $make_status
                echo
                echo "ERROR: Compilation of mega2 failed.  Get help"
                exit $make_status
            fi
            mega2prog=srcdir/mega2_${OSTYPE}
            cd ..
            if [[ ${SAVE:-""} != "" ]]; then
                echo
                echo saving $mega2prog to mega2_bin/$pgm
                rm -f mega2_bin/$pgm
                cp -p $mega2prog mega2_bin/$pgm
            fi
        else 
            echo "ERROR: No pre-built Makedefs_${OSTYPE}.template file"
            echo "To compile mega2, "
            echo "   Please create a Makedefs_${OSTYPE}.template file inside the srcdir folder"
            echo "   and run install.sh again."
            echo "   You may using an existing Makedefs template file as an example."
            exit
        fi
    fi
done

################################################################
#                      INSTALL
################################################################

verify

if [ ! -d $mega2_path ]; then
    echo "Directory $mega2_path does not exist."
    echo "Please run install.sh again"
    exit
fi


for scr in $scripts; do
    if [[ $scr == "mega2" || $scr == "mega2compile" ]]; then

	if [ -x $mega2prog ]; then
	    echo "Installing  $mega2prog for $OSTYPE version $v"
	fi

        if [[ $MEGA2_OVERWRITE_BINARIES ]] ; then
            /bin/cp -p $mega2prog $mega2_path"/mega2_${VERSION}_${OSTYPE}"
	elif [[ -e $mega2_path"/mega2_${VERSION}_${OSTYPE}" ||  -h $mega2_path"/mega2_${VERSION}_${OSTYPE}" ]]; then
	    echo "A file named mega2_${VERSION}_${OSTYPE} already exists in $mega2_path"
	    echo -n "Remove this file [y/n]? " 
            read answer
	    if [[ $answer == 'y' ]]; then
		/bin/rm -f $mega2_path"/mega2_${VERSION}_${OSTYPE}"
		/bin/cp -p $mega2prog $mega2_path"/mega2_${VERSION}_${OSTYPE}"
	    else
		echo "Not removing $mega2_path/mega2_${VERSION}_${OSTYPE}"
	    fi
        else
		/bin/cp -p $mega2prog $mega2_path"/mega2_${VERSION}_${OSTYPE}"
	fi

        if [[ $MEGA2_OVERWRITE_BINARIES ]] ; then
            /bin/cp -p $mega2prog $mega2_path"/mega2"
	elif [[ -e $mega2_path"/mega2" || -h $mega2_path"/mega2" ]]; then
	    echo "A file named mega2 already exists in $mega2_path"
	    echo -n "Remove this file [y/n]? "
            read answer
	    if [[ $answer == 'y' ]]; then
		/bin/rm -f $mega2_path"/mega2"
		/bin/cp -p $mega2prog $mega2_path"/mega2"
	    else
		echo "Not removing $mega2_path/mega2_${VERSION}_${OSTYPE}"
	    fi
        else
		/bin/cp -p $mega2prog $mega2_path"/mega2"
	fi


    elif [[ $scr == "scripts" ]]; then

	if [[ ! $perl_path || $perl_path == "" ]]; then
	    echo "ERROR: Could not find path where perl is installed."
	    echo "Failure installing perl scripts."
	else 
	    for pl in $perl_scripts; do
		if [ ! -f $pl.src ]; then
		    echo "WARNING: Source file $pl.src for perl script $pl not found,"
		    echo Skipping ..
		else
		    if [ -f "$mega2_path/$pl" ]; then 
			/bin/rm "$mega2_path/$pl";
		    fi
		
		    echo "#!"$perl_path" -w" > $mega2_path"/"$pl
		    datestr=`date`;
		    echo "#Installed: $datestr" >> $mega2_path"/"$pl
		    echo "#" >> $mega2_path"/"$pl
		    echo "#=====end automatic header=======" >> $mega2_path"/"$pl
		    echo "#" >> $mega2_path"/"$pl

		    cat $pl".src" >> $mega2_path"/"$pl
		    chmod +x $mega2_path"/"$pl

                    if [ ! -x $mega2_path"/"$pl ]; then
                        echo "Failed to install perl script $pl in $mega2_path"
                    else 
                        echo "Installed perl script $pl in $mega2_path"
                    fi
		fi
	    done
	fi
	
	if [[ ! $python_path || $python_path == "" ]]; then
	    echo "ERROR: Could not find path where python is installed."
	    echo "Failure installing python scripts"
	else
	    for py in $python_scripts; do
		if [ ! -f $pl.src ]; then
		    echo "WARNING: Source file $pl.src for perl script $pl not found,"
		    echo Skipping ..
		else
		    if [ -f "$mega2_path/$py" ]; then 
			/bin/rm "$mega2_path/$py";
		    fi
		
		    echo "#!"$python_path > $mega2_path"/"$py
		    datestr=`date`;
		    echo "#Installed: $datestr" >> $mega2_path"/"$py
		    echo "#" >> $mega2_path"/"$py
		    echo "#=====end automatic header=======" >> $mega2_path"/"$py
		    echo "#" >> $mega2_path"/"$py
		    cat $py".src" >> $mega2_path"/"$py
		    chmod +x $mega2_path"/"$py

                    if  [ ! -x $mega2_path"/"$py ]; then
                        echo "Failed to install $py in $mega2_path"
                    else 
                        echo "Installed $py in $mega2_path"
                    fi	    
		fi	
	    done
	fi
    fi
done

echo
echo "Please remember to also install the 'nplplot' and 'genetics' R packages into R if necessary"
echo "Please consult the Mega2 documentation for details"
