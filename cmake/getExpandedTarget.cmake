function(expand_target_all_targets target_names expanded_list)
    foreach(target_name ${target_names})
        if(TARGET ${target_name} AND NOT ${target_name} IN_LIST expanded_list)
            list(APPEND target_names_expanded ${target_name})
            unset(interface_libs)
            unset(private_libs)
            unset(imported_lib)
            unset(lib_type)
            unset(lib_impt)
            get_target_property(lib_type ${target_name} TYPE)
            get_target_property(lib_impt ${target_name} IMPORTED)
            if(lib_impt)
                if(NOT lib_type MATCHES "INTERFACE" OR CMAKE_VERSION VERSION_GREATER_EQUAL 3.19)
                    # The location property can only be read on imported targets, or interface imported when cmake > 3.19
                    get_target_property(imported_lib ${target_name} LOCATION)
                endif()
            endif()
            if(NOT lib_type MATCHES "INTERFACE")
                get_target_property(private_libs ${target_name} LINK_LIBRARIES)
            endif()
            get_target_property(interface_libs ${target_name} INTERFACE_LINK_LIBRARIES)
            list(FILTER imported_lib EXCLUDE REGEX "NOTFOUND|(-o)$")
            list(FILTER private_libs EXCLUDE REGEX "NOTFOUND|(-o)$")
            list(FILTER interface_libs EXCLUDE REGEX "NOTFOUND|(-o)$")
            foreach(elem ${imported_lib};${private_libs};${interface_libs})
                string(REGEX REPLACE "([\$<]+[A-Za-z]+:[A-Za-z]+[:>]+)|:>|>" "" elem_stripped "${elem}")
                if(NOT TARGET ${elem_stripped})
                    continue()
                endif()
                if(${elem_stripped} IN_LIST target_names)
                    continue()
                endif()
                if(${elem} IN_LIST expanded_list)
                    continue()
                endif()
                if(${elem_stripped} IN_LIST target_names_expanded)
                    continue()
                endif()
                unset(recursed_list) # Otherwise this one grows for each elem
                expand_target_all_targets(${elem_stripped} recursed_list)
                list(REMOVE_DUPLICATES "recursed_list")
                foreach(rec ${recursed_list})
                    if(NOT TARGET ${rec})
                        continue()
                    endif()
                    if(${rec} IN_LIST target_names)
                        continue()
                    endif()
                    if(${rec} IN_LIST expanded_list)
                        continue()
                    endif()
                    if(${rec} IN_LIST target_names_expanded)
                        continue()
                    endif()
                    list(APPEND target_names_expanded ${rec})
                endforeach()
            endforeach()
        endif()
    endforeach()

    # Remove duplicates in a way that retains linking order, i.e. keep last occurrence
    if(target_names_expanded)
        list(REVERSE target_names_expanded)
        list(REMOVE_DUPLICATES "target_names_expanded")
        list(REVERSE "target_names_expanded")
        list(APPEND ${target_names_expanded} "${${expanded_list}}")
    endif()
    list(APPEND ${expanded_list} "${target_names_expanded}")
    set(${expanded_list} "${${expanded_list}}" PARENT_SCOPE)
endfunction()
