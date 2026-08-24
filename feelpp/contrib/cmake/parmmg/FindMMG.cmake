# ParMmg 1.5.0 searches for a pre-built MMG library when DOWNLOAD_MMG is
# disabled. Feel++ builds MMG immediately before ParMmg, so expose that target
# through the variables expected by ParMmg's existing integration code.
if ( NOT TARGET Mmg::libmmg_so )
  set( MMG_FOUND OFF )
  return()
endif()

set( MMG_FOUND ON )
set( MMG_LIBRARIES Mmg::libmmg_so )
set( MMG_INCLUDE_DIRS "${MMG_BINARY_DIR}/include" )
set( MMG_LIBRARY_DIRS "" )
