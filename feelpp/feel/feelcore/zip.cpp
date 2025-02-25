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
        //throw std::runtime_error( fmt::format("Could not open zip file {}", zipFilePath) );
        LOG(ERROR) << fmt::format( "Could not open zip file {}", zipFilePath ) << std::endl;
        return false;
    }

    int numEntries = zip_get_num_entries( archive, 0 );
    for ( int i = 0; i < numEntries; ++i )
    {
        const char* entryName = zip_get_name( archive, i, 0 );
        if ( !entryName )
        {
            zip_close( archive );
            LOG(ERROR) << fmt::format( "Could not get entry name " ) << std::endl;
            return false;
        }

        fs::path extractionPath = fs::path(extractionDir) / fs::path( entryName );
        VLOG(2) << fmt::format( "Extracting {} to {}", entryName, extractionPath.string() ) << std::endl;
        VLOG(2) << fmt::format( "has {} a filename : {}", entryName, fs::path( entryName ).has_filename() ) << std::endl;
        if ( !fs::path( entryName ).has_filename() || fs::path( entryName ) == fs::path(".") )
        {
            VLOG(2) << fmt::format( "Creating directory {}", extractionPath.string() ) << std::endl;
            fs::create_directories( extractionPath );
            CHECK( fs::is_directory( extractionPath ) ) << fmt::format( "Could not create directory {}", extractionPath.string() );
            continue;
        }

        zip_file* file = zip_fopen_index( archive, i, 0 );
        if ( file && fs::is_directory( extractionPath ) )
        {
            zip_fclose( file );
            continue;
        }
        else if ( !file )
        {
            zip_close( archive );
            LOG(ERROR) << fmt::format( "Could not open file {} in zip file", entryName ) << std::endl;
            return false;
        }
        else
        {
            FILE* outFile = fopen( extractionPath.string().c_str(), "wb" );
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
                fwrite( buf, 1, bytesRead, outFile );
            }

            fclose( outFile );
            zip_fclose( file );
        }
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
