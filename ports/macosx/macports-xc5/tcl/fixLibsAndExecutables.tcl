#!/usr/bin/tclsh

# does a depth first search in the directory dirs 
# to find executables and dynamic libraries
proc findLibsAndExecutables {{dirs .}} {
    set lst [list]
    while {[llength $dirs] != 0} {
       set name [lindex $dirs 0]
       set dirs [concat [glob -nocomplain -directory [lindex $dirs 0] *] [lrange $dirs 1 end]]
       #puts $name

       if { [file extension $name] == ".dylib" 
       || [file isfile $name ] && [file executable $name]
       } then {
           lappend lst $name
       }
    }

    return $lst
}

# Look for the file fname
# in the directory tree dirs 
proc findRelativePath {dirs fname} {
    while {[llength $dirs] != 0} {
       set name [lindex $dirs 0]
       set dirs [concat [glob -nocomplain -directory [lindex $dirs 0] *] [lrange $dirs 1 end]]
       #puts $name

       if { [file tail $name] == $fname } then {
           return $name
       }
    }

    return ""
}

# Tries to fix bad library references in libraries and executables of a package built from sources
# This procedure is to be used in the post-destroot step
proc fixLibsAndExecutables {workpath destroot} {
    # Find libraries and executables of the current projet
    set files [findLibsAndExecutables $destroot]
    #puts $files

    # iterate over the libraries and executables
    foreach fname $files {
        puts "File: $fname"
        set command "otool -L $fname"
        if {[catch {eval exec $command} result] == 0} { 
            #puts $result
            set libs [split $result "\n"]
            #puts $libs
    
            # iterate over the libraries referenced by the current file 
            foreach l $libs {
                #puts "$l is item number in list x"
                # trim spaces
                set ltr [string trim $l]

                # if we find a reference to a library that has not been changed
                # or a reference to a library that is not an absolute path
                if { [string first $workpath $ltr] == 0 || [string first "/" $ltr] != 0 } {
                    # extract the library name
                    set spl [split $ltr " "]
	            set lchange [lindex $spl 0]
                    #puts $lchange 
                    #puts [file tail $lchange]

                    # find the correct path to the library
                    set lnew [findRelativePath $destroot [file tail $lchange]]
                    if { $lnew != "" } {
                        #puts $lnew

                        set lg [string length $destroot]
                        set lnew [string replace $lnew 0 [expr "$lg - 1"]]

                        # Fix the path
                        puts "Fixing entry: $lchange -> $lnew"
                        set command "install_name_tool -change $lchange $lnew $fname"
                        eval exec $command
                    }
                }
            }
        } else { 
            puts "Error while using libtool"
        } 
    }
}

#set destroot "/opt/local/var/macports/software/petsc/opt/local/lib/"
#set destroot "/opt/local/var/macports/build/_Users_aancel_code_feelpp_ports_macosx_macports-xc5_science_feelpp/feel++/work/destroot"
#set workpath "/opt/local/var/macports/build/_Users_aancel_code_feelpp_ports_macosx_macports-xc5_science_petsc/petsc/work/"
#set workpath "/opt/local/var/macports/build/_Users_aancel_code_feelpp_ports_macosx_macports-xc5_science_feelpp/feel++/work"

#fixLibsAndExecutables $workpath $destroot

