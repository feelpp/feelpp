# Microsoft Developer Studio Project File - Name="Library" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=Library - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "Library.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "Library.mak" CFG="Library - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "Library - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "Library - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "Library - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /I "c:/spirit161" /I "c:/boost_1_31_0" /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /Zm800 /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "Library - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GR /GX /Zi /Od /I "c:\spirit161" /I "..\..\.." /I "c:/spirit161" /I "c:/boost_1_31_0" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /FD /GZ /Zm800 /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ENDIF 

# Begin Target

# Name "Library - Win32 Release"
# Name "Library - Win32 Debug"
# Begin Group "Achive Implementations"

# PROP Default_Filter "hpp"
# Begin Source File

SOURCE=..\..\..\boost\archive\archive_exception.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\basic_archive.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\binary_iarchive.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\binary_oarchive.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\extended_type_info.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\polymorphic_binary_iarchive.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\polymorphic_binary_oarchive.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\polymorphic_iarchive.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\polymorphic_oarchive.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\polymorphic_text_iarchive.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\polymorphic_text_oarchive.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\polymorphic_xml_iarchive.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\polymorphic_xml_oarchive.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\text_iarchive.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\text_oarchive.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\xml_iarchive.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\xml_oarchive.hpp
# End Source File
# End Group
# Begin Group "Serialization User Headers"

# PROP Default_Filter "hpp"
# Begin Source File

SOURCE=..\..\..\boost\serialization\access.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\serialization\base_object.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\serialization\binary_object.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\serialization\export.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\serialization\extended_type_info.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\serialization\extended_type_info_no_rtti.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\serialization\extended_type_info_typeid.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\serialization\is_abstract.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\serialization\level.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\serialization\level_enum.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\serialization\split_free.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\serialization\split_member.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\serialization\tracking.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\serialization\tracking_enum.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\serialization\traits.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\serialization\type_info_implementation.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\serialization\version.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\serialization\void_cast.hpp
# End Source File
# End Group
# Begin Group "Serialization Implementations"

# PROP Default_Filter "hpp"
# Begin Source File

SOURCE=..\..\..\boost\serialization\collection_traits.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\serialization\collections_load_imp.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\serialization\collections_save_imp.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\serialization\deque.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\serialization\hash_map.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\serialization\hash_set.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\serialization\list.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\serialization\map.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\serialization\nvp.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\serialization\optional.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\serialization\scoped_ptr.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\serialization\set.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\serialization\shared_count.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\serialization\shared_ptr.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\serialization\slist.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\serialization\string.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\serialization\utility.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\serialization\vector.hpp
# End Source File
# End Group
# Begin Group "Archive Developer Headers"

# PROP Default_Filter "hpp"
# Begin Source File

SOURCE=..\..\..\boost\archive\basic_binary_iarchive.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\basic_binary_iprimitive.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\basic_binary_oarchive.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\basic_binary_oprimitive.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\basic_text_iarchive.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\basic_text_iprimitive.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\basic_text_oarchive.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\basic_text_oprimitive.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\basic_xml_archive.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\basic_xml_iarchive.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\basic_xml_oarchive.hpp
# End Source File
# End Group
# Begin Group "Archive Detail"

# PROP Default_Filter "hpp"
# Begin Source File

SOURCE=..\..\..\boost\archive\detail\archive_pointer_iserializer.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\detail\archive_pointer_oserializer.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\detail\basic_archive_pointer_iserializer.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\detail\basic_archive_pointer_oserializer.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\detail\basic_iarchive.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\detail\basic_iserializer.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\detail\basic_oarchive.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\detail\basic_oserializer.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\detail\basic_pointer_iserializer.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\detail\basic_pointer_oserializer.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\detail\basic_serializer.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\detail\basic_serializer_map.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\detail\common_archive.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\detail\common_iarchive.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\detail\common_oarchive.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\detail\interface_iarchive.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\detail\interface_oarchive.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\detail\iserializer.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\detail\known_archive_types.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\detail\known_archive_types_fwd.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\detail\oserializer.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\detail\polymorphic_iarchive_impl.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\detail\polymorphic_oarchive_impl.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\serialization\serialization.hpp
# End Source File
# End Group
# Begin Group "Utility Headers"

# PROP Default_Filter "hpp"
# Begin Source File

SOURCE=..\..\..\boost\archive\add_facet.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\codecvt_null.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\dinkumware.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\serialization\force_include.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\pfto.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\smart_cast.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\state_saver.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\static_warning.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\strong_typedef.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\utf8_codecvt_facet.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\wcslen.hpp
# End Source File
# End Group
# Begin Group "Templates"

# PROP Default_Filter "ipp"
# Begin Source File

SOURCE=..\..\..\boost\archive\impl\archive_pointer_iserializer.ipp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\impl\archive_pointer_oserializer.ipp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\impl\basic_binary_iprimitive.ipp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\impl\basic_binary_iprimitive.ipp.bak
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\impl\basic_binary_oprimitive.ipp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\impl\basic_text_iarchive.ipp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\impl\basic_text_iprimitive.ipp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\impl\basic_text_oprimitive.ipp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\impl\basic_xml_iarchive.ipp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\impl\basic_xml_oarchive.ipp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\impl\text_iarchive_impl.ipp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\impl\text_oarchive_impl.ipp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\impl\text_wiarchive_impl.ipp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\impl\text_woarchive_impl.ipp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\impl\xml_iarchive_impl.ipp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\impl\xml_oarchive_impl.ipp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\impl\xml_wiarchive_impl.ipp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\impl\xml_woarchive_impl.ipp
# End Source File
# End Group
# Begin Group "Iterators"

# PROP Default_Filter "hpp"
# Begin Source File

SOURCE=..\..\..\boost\archive\iterators\base64_exception.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\iterators\base64_from_binary.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\iterators\binary_from_base64.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\iterators\escape.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\iterators\insert_linebreaks.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\iterators\istream_iterator.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\iterators\mb_from_wchar.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\iterators\ostream_iterator.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\iterators\remove_whitespace.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\iterators\transform_width.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\iterators\unescape.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\iterators\wchar_from_mb.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\iterators\xml_escape.hpp
# End Source File
# Begin Source File

SOURCE=..\..\..\boost\archive\iterators\xml_unescape.hpp
# End Source File
# End Group
# Begin Group "Source"

# PROP Default_Filter "cpp"
# Begin Source File

SOURCE=..\src\basic_archive.cpp
# End Source File
# Begin Source File

SOURCE=..\src\basic_iarchive.cpp
# End Source File
# Begin Source File

SOURCE=..\src\basic_oarchive.cpp
# End Source File
# Begin Source File

SOURCE=..\src\basic_serializer_map.cpp
# End Source File
# Begin Source File

SOURCE=..\src\basic_text_iprimitive.cpp
# End Source File
# Begin Source File

SOURCE=..\src\basic_text_oprimitive.cpp
# End Source File
# Begin Source File

SOURCE=..\src\basic_xml_archive.cpp
# End Source File
# Begin Source File

SOURCE=..\src\basic_xml_grammar.ipp
# End Source File
# Begin Source File

SOURCE=..\src\binary_iarchive.cpp
# End Source File
# Begin Source File

SOURCE=..\src\binary_oarchive.cpp
# End Source File
# Begin Source File

SOURCE=..\src\codecvt_null.cpp
# End Source File
# Begin Source File

SOURCE=..\src\extended_type_info.cpp
# End Source File
# Begin Source File

SOURCE=..\src\extended_type_info_no_rtti.cpp
# End Source File
# Begin Source File

SOURCE=..\src\extended_type_info_typeid.cpp
# End Source File
# Begin Source File

SOURCE=..\src\polymorphic_iarchive.cpp
# End Source File
# Begin Source File

SOURCE=..\src\polymorphic_oarchive.cpp
# End Source File
# Begin Source File

SOURCE=..\src\text_iarchive.cpp
# End Source File
# Begin Source File

SOURCE=..\src\text_oarchive.cpp
# End Source File
# Begin Source File

SOURCE=..\src\void_cast.cpp
# End Source File
# Begin Source File

SOURCE=..\src\xml_grammar.cpp
# End Source File
# Begin Source File

SOURCE=..\src\xml_iarchive.cpp
# End Source File
# Begin Source File

SOURCE=..\src\xml_oarchive.cpp
# End Source File
# End Group
# End Target
# End Project
