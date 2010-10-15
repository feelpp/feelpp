file(GLOB Eigen_directory_files "*")
foreach(f ${Eigen_directory_files})
  if( NOT f MATCHES ".txt" AND NOT IS_DIRECTORY ${f} )
    list(APPEND Eigen_directory_files_to_install ${f})
  endif()
endforeach(f ${Eigen_directory_files})

install(FILES
  ${Eigen_directory_files_to_install}
  DESTINATION include/Eigen COMPONENT Devel
  )

add_subdirectory(src)
