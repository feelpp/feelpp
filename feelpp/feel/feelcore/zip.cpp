#include <cstdio>
#include <cstdlib>
#include <string>
#include <zip.h>
#include <fmt/core.h>
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/zip.hpp>

namespace Feel
{
/**
 * @brief Extract a zip file to a directory
 *
 * @param zipFilePath path to the zip file
 * @param extractionDir directory where to extract the zip file
 * @return true if the extraction was successful
 * @return false otherwise
 */
bool extractZipFile( const std::string& zipFilePath, const std::string& extractionDir )
{
    LOG(INFO) << fmt::format( "Extracting {} to {}", zipFilePath, extractionDir ) << std::endl;
    zip* archive = zip_open( zipFilePath.c_str(), 0, nullptr );
    if ( !archive )
    {
        LOG(ERROR) << fmt::format( "Could not open zip file {}", zipFilePath ) << std::endl;
        return false;
    }

    int numEntries = zip_get_num_entries( archive, 0 );
    // Browse all files in the archive
    for ( int i = 0; i < numEntries; ++i )
    {
        // Get filename at position i in the archive
        const char* entryName = zip_get_name( archive, i, 0 );
        if ( !entryName )
        {
            zip_close( archive );
            LOG(ERROR) << fmt::format( "Could not get entry name " ) << std::endl;
            return false;
        }

        fs::path extractionPath = fs::path(extractionDir) / fs::path( entryName );
        std::cout << fmt::format( "Extracting {} to {}", entryName, extractionPath.string() ) << std::endl;

        // If it's a folder
        struct zip_stat st;
        zip_stat_index(archive, i, 0, &st);
        bool isDirectory = (st.name[strlen(st.name) - 1] == '/');

        if (isDirectory)
        {
            VLOG(2) << fmt::format( "Creating directory {}", extractionPath.string() ) << std::endl;
            fs::create_directories( extractionPath );
            continue;
        }

        // If it's a file
        zip_file* file = zip_fopen_index( archive, i, 0 );
        if ( !file )
        {
            zip_close( archive );
            LOG(ERROR) << fmt::format( "Could not open file {} in zip file", entryName ) << std::endl;
            return false;
        }

        fs::create_directories( extractionPath.parent_path() ); // create extractionPath's parent folder if it doesn't exist
        std::ofstream outFile(extractionPath, std::ios::binary);
        if ( !outFile )
        {
            zip_fclose( file );
            zip_close( archive );
            LOG(ERROR) << fmt::format( "Could not open output file {}", extractionPath.string() ) << std::endl;
            return false;
        }

        zip_int64_t bytesRead;
        char buf[8192];
        while ( ( bytesRead = zip_fread( file, buf, sizeof( buf ) ) ) > 0 )
        {
            outFile.write(buf, bytesRead);
        }

        outFile.close();
        zip_fclose( file );
    }

    zip_close( archive );
    return true;
}

/**
 * @brief Remove all files and directories in a directory
 *
 * @param extractionDir directory to clean up
 */
void cleanupTemporaryDirectory( const std::string& extractionDir )
{
    #if 0
    // Iterate through the files and directories in the extraction directory
    for ( const auto& entry : fs::directory_iterator( extractionDir ) )
    {
        try
        {
            if ( fs::is_directory( entry ) )
            {
                fs::remove_all( entry ); // Remove directories recursively
            }
            else
            {
                fs::remove( entry ); // Remove individual files
            }
        }
        catch ( const std::exception& e )
        {
            std::cerr << "Error while cleaning up: " << e.what() << std::endl;
        }
    }
    #endif
    fs::remove_all( extractionDir );
}
} // namespace Feel
